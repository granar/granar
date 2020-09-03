
#' @title GRANAR: Generate root anatomy in R
#'
#' This function creates a 2D root cross section anatomy based on anatomical features translated into the parameter of the model. It is the core of GRANAR.
#' @param path Path to the XML file containing the different parameters for the simulation. Not needed if 'parameter' is set.  Default = NULL.
#' @param parameters Table with the different parameters. Not needed if 'path' is set.  Default = NULL.
#' @param verbatim Display textual information aboutt he simulation. Default = NULL.
#' @param maturity_x if True Metaxylem vessels are tagged as stele cells. Default = FALSE.
#' @keywords root anatomy
#' @author Adrien Heymans and Guillaume Lobet
#' @export
#' @import tidyverse
#' @import plyr
#' @import deldir
#' @import alphahull
#' @import sp
#' @examples
#'
#' root_cross_section <- create_anatomy(path = "PATH_TO_XLM_FILE")
#'
#' # OR
#'
#' params <- read_param_xml(path = "PATH_TO_XLM_FILE")
#' root_cross_section <- create_anatomy(parameters = params)
#'
#' # To visualize the result
#'
#' plot_anatomy(root_cross_section)

create_anatomy <- function(path = NULL,  # PAth
                             parameters = NULL,
                             verbatim = F,
                             maturity_x = F){

  # Return NULL is no parameters are specified
  if( is.null(path) & is.null(parameters)){
    warning("Please specify a parameter set for the simulation")
    return(NULL)
  }
  if(!is.null(path)){
    params <- read_param_xml(path)
  }
  if(!(is.null(parameters))){
    params <- parameters
    if(nrow(params[params$name == "planttype",]) == 0){
      warning(paste0("Could not find the 'planttype' information in the parameter input"))
      return(NULL)
    }
    # Quality control
    to_find <- c("planttype", "randomness", "xylem", "phloem", "stele", "endodermis", "exodermis", "epidermis", "aerenchyma", "pericycle", "cortex", "n_passage_cell")
    for(tf in to_find){
      if (nrow(params[params$name == tf,]) == 0){
        warning(paste0("Could not find the '",tf,"' information in the parameter input"))
      }
    }
    cols_to_find <- c("name", "type", "value")
    for(ctf in cols_to_find){
      if (is.null(params[[ctf]])){
        warning(paste0("Could not find the '",ctf,"' column in the parameter input"))
        return(NULL)
      }
    }

  }

  t_1 <- Sys.time()
  t1 <- proc.time()

  # PARAMETERS -----
  if(verbatim) message("Loading parameters")
  plant_type <- params$value[params$name == "planttype"]
  random_fact <- params$value[params$name == "randomness"] / 10 * params$value[params$name == "stele" & params$type == "cell_diameter"]
  stele_diameter <- params$value[params$name == "stele" & params$type == "layer_diameter"]
  n_aerenchyma_files <- params$value[params$name == "aerenchyma" & params$type == "n_files"]
  proportion_aerenchyma <- params$value[params$name == "aerenchyma" & params$type == "proportion"]
  n_xylem_files <- params$value[params$name == "xylem" & params$type == "n_files"]
  proto_meta_ratio <- params$value[params$name == "xylem" & params$type == "ratio"]
  n_proto_xylem <- round(n_xylem_files*proto_meta_ratio)
  if(length(params$value[params$name == "n_passage_cell"]) == 0L){
    n_passar <- 0
  }else{n_passar <- params$value[params$name == "n_passage_cell"]}

  t2 <- proc.time()

  # params$value[params$name == "stele" & params$type == "cell_diameter"] <- params$value[params$name == "stele" & params$type == "cell_diameter"]-params$value[params$name == "stele" & params$type == "cell_diameter"]*0.084
  # params$value[params$name == "pericycle" & params$type == "cell_diameter"] <- params$value[params$name == "pericycle" & params$type == "cell_diameter"]-params$value[params$name == "pericycle" & params$type == "cell_diameter"]*0.17
  # params$value[params$name == "endodermis" & params$type == "cell_diameter"] <- params$value[params$name == "endodermis" & params$type == "cell_diameter"]-params$value[params$name == "endodermis" & params$type == "cell_diameter"]*0.449
  # params$value[params$name == "cortex" & params$type == "cell_diameter"] <- params$value[params$name == "cortex" & params$type == "cell_diameter"]-params$value[params$name == "cortex" & params$type == "cell_diameter"]*0.105
  # params$value[params$name == "exodermis" & params$type == "cell_diameter"] <- params$value[params$name == "exodermis" & params$type == "cell_diameter"]- params$value[params$name == "exodermis" & params$type == "cell_diameter"]*0.0916
  # params$value[params$name == "epidermis" & params$type == "cell_diameter"] <- params$value[params$name == "epidermis" & params$type == "cell_diameter"]-params$value[params$name == "epidermis" & params$type == "cell_diameter"]*0.225

  #//////////////////////////////////////////////////////////////////////////////////////////////////
  # INITIALIZE LAYERS -----

  if(verbatim) message("Creating cell layers")

  layers <- params %>%
    filter(type %in% c("cell_diameter","n_layers","order")) %>%
    spread(type, value) %>%
    filter(!is.na(n_layers)) %>%
    arrange(order)

  # Create and "outside" layer to serve as boundary for the voronoi algorithm.
  layers <- rbind(layers, data.frame(name="outside",
                                     n_layers=2,
                                     cell_diameter=layers$cell_diameter[layers$name == "epidermis"]* 1,
                                     order = max(layers$order)+1))

  # Get the number of cell layers for the stele
  layers$n_layers[layers$name == "stele"] <- round((stele_diameter/2) / layers$cell_diameter[layers$name == "stele"]) #
  #layers$size[layers$name == "stele"] <- diam_stele


  # Get one row per actual cell layer
  all_layers <- NULL
  for(i in c(1:nrow(layers))){
    for(j in c(1:layers$n_layers[i])){
      all_layers <- rbind(all_layers, layers[i,])
    }
  }

  all_layers$radius <- all_layers$cell_diameter / 2
  all_layers$perim <- all_layers$radius * 2 * pi
  all_layers$n_cell <- 1
  all_layers$angle_inc <- 0

  all_layers$radius[1] <- 0
  multi <- 1

  for(i in c(2:nrow(all_layers))){
    # Update radius
    all_layers$radius[i] <- all_layers$radius[i-1] + all_layers$cell_diameter[i-1] / 2 + all_layers$cell_diameter[i] / 2

    if(all_layers$name[i] == "pericyle" & multi == 1){
      all_layers$radius[i] <- stele_diameter/2 + all_layers$cell_diameter[i] / 2
      multi <- 2 # in case there is more than one pericycle layer
    }
    # Update perimeter
    all_layers$perim[i] <- all_layers$radius[i] * 2 * pi

    # Update number of cells in the layers
    all_layers$n_cell[i] <- round(all_layers$perim[i] / all_layers$cell_diameter[i])

    # Update the mean angle between cells
    all_layers$angle_inc[i] <- 2 * pi / all_layers$n_cell[i]
  }


  t3 <- proc.time()



  #//////////////////////////////////////////////////////////////////////////////////////////////////
  # CREATE CELLS ------

  if(verbatim) message("Creating cells")

  center <- max(all_layers$radius)
  all_cells <- NULL
  k <- 1
  for(i in c(1:nrow(all_layers))){
    radius <- all_layers$radius[i]
    if(all_layers$angle_inc[i] > 0){
      angles <- seq(from = 0, to = (2*pi), by = all_layers$angle_inc[i])[-1]
    }else{
      angles <- 0
    }
    k1 <- k+all_layers$n_cell[i]-1
    ks <- c(k:k1)
    k <- k1+1

    if(all_layers$name[i] == "outside"){
      x <- center + (radius * cos(angles))
      y <- center + (radius * sin(angles))
    }else if(all_layers$name[i] == "stele"){
      x <- center + (radius * cos(angles)) + runif(all_layers$n_cell[i], -random_fact, random_fact)#* runif(all_layers$n_cell[i], 1-(random_fact*2), 1+(random_fact*2))
      y <- center + (radius * sin(angles)) + runif(all_layers$n_cell[i], -random_fact, random_fact)
    }else if(all_layers$name[i] == "cortex"){
      x <- center + (radius * cos(angles)) + runif(all_layers$n_cell[i], -random_fact*3, random_fact*3)#* runif(all_layers$n_cell[i], 1-(random_fact*2), 1+(random_fact*2))
      y <- center + (radius * sin(angles)) + runif(all_layers$n_cell[i], -random_fact*3, random_fact*3)##* runif(all_layers$n_cell[i], 1-(random_fact*2), 1+(random_fact*2))
    }else{
      x <- center + (radius * cos(angles)) + runif(all_layers$n_cell[i], -random_fact, random_fact)#* runif(all_layers$n_cell[i], 1-random_fact, 1+random_fact)
      y <- center + (radius * sin(angles)) + runif(all_layers$n_cell[i], -random_fact, random_fact)#* runif(all_layers$n_cell[i], 1-random_fact, 1+random_fact)
    }

    all_cells <- rbind(all_cells, data.frame(
      angle = angles,
      radius = radius,
      x = x,
      y = y,
      id_layer = i,
      id_cell = ks,
      type = all_layers$name[i],
      order = all_layers$order[i]
    )
    )
  }

  # all_cells <- all_cells_bis
  summary_cells <- ddply(all_cells, .(type), summarise, n_cells = length(angle))

  t4 <- proc.time()

  #//////////////////////////////////////////////////////////////////////////////////////////////////
  # CREATE XYLEM VESSELS -----
  # Create the xylem files
  # Get the extremes

  if(verbatim) message("Creating xylem and phloem vessels")

  all_cells$id_group <- 0

  # Create vascular vessels for dicot or monocot
  if(plant_type == 2){ # DICOT
    xyl <- data.frame(r=numeric(2), d=numeric(2))
    xyl$r <- c(0, max(all_cells$radius[all_cells$type == "stele"]))
    xyl$d <- c(params$value[params$type == "max_size" & params$name == "xylem"], layers$cell_diameter[layers$name == "stele"])

    # Get the cells in between
    fit <- lm(d ~ r, data=xyl)$coefficients
    rnew <- xyl$r[1]
    i <- 1
    rmax <- xyl$r[2]
    dmin <- xyl$d[2]
    keep_going <- T
    while(keep_going){
      xyl <- xyl %>% arrange(r)

      rnew <- xyl$r[i] + xyl$d[i] #+ (xyl$r[2]/10)
      dnew <- fit[1] + rnew*fit[2]
      while(rnew+(dnew/2) > rmax-(dmin/2)){
        rnew <- rnew - 0.05
        dnew <- dnew - 0.05
        keep_going = F
      }
      xyl <- rbind(xyl, data.frame(r = rnew,d = dnew))
      i <- i+1
    }
    xyl <- xyl %>% arrange(r)%>%
      filter(d > 0)
    while(xyl$d[nrow(xyl)] >= xyl$d[nrow(xyl)-1]){
      xyl$d[nrow(xyl)] <- xyl$d[nrow(xyl)] - 0.04
      xyl$r[nrow(xyl)] <- xyl$r[nrow(xyl)] - 0.02
    }
    all_xylem <- NULL
    i <- 1
    angle_seq <- seq(from = 0, to = (2*pi), by = (2 * pi) / n_xylem_files)
    x <- center + (xyl$r[1] * cos(angle_seq[1]))
    y <- center + (xyl$r[1] * sin(angle_seq[1]))
    all_xylem <- rbind(all_xylem, data.frame(x = x,
                                             y = y,
                                             d = xyl$d[1],
                                             angle = angle_seq[1],
                                             id_group = i))
    all_cells <- rbind(all_cells, data.frame(
      angle = angle_seq[1],
      radius = xyl$r[1],
      x = x,
      y = y,
      id_layer = 20,
      id_cell = 1,
      type = "xylem",
      order = 1.5,
      id_group = i
    ))

    i <- i+1
    for(angle in angle_seq){
      x <- center + (xyl$r[-1] * cos(angle))
      y <- center + (xyl$r[-1] * sin(angle))
      all_xylem <- rbind(all_xylem, data.frame(x = x,
                                               y = y,
                                               d = xyl$d[-1],
                                               angle = angle,
                                               id_group = i))
      all_cells <- rbind(all_cells, data.frame(
        angle = angle,
        radius = xyl$r[-1],
        x = x,
        y = y,
        id_layer = 20,
        id_cell = 1,
        type = "xylem",
        order = 1.5,
        id_group = i
      )
      )
      i <- i+1
    }

  }else if(plant_type == 1){ # MONOCOT
    if(n_xylem_files == 1){ # One metaxylem in the center of the stele
      r <- 0
      xyl <- data.frame(r = r,
                        d = params$value[params$type == "max_size" & params$name == "xylem"])
    }else{

      #modification 05/01
      r= max(all_cells$radius[all_cells$type == "stele"]) - (params$value[params$type == "cell_diameter" & params$name == "stele"])*1.5 -
        (params$value[params$type == "max_size" & params$name == "xylem"])/2
      xyl <- data.frame(r = r,
                        d = params$value[params$type == "max_size" & params$name == "xylem"])
    }
    #  xyl <- data.frame(r = max(all_cells$radius[all_cells$type == "stele"]) - (params$value[params$type == "max_size" & params$name == "xylem"])/2,
    #                  d = params$value[params$type == "max_size" & params$name == "xylem"])
    all_xylem <- NULL
    angle_seq <- seq(from = 0, to = (2*pi), by = (2 * pi) / n_xylem_files)
    i <- 1
    for(angle in angle_seq){
      x <- center + (xyl$r[1] * cos(angle))
      y <- center + (xyl$r[1] * sin(angle))
      all_xylem <- rbind(all_xylem, data.frame(x = x,
                                               y = y,
                                               d = xyl$d[1],
                                               angle = angle,
                                               id_group = i))
      all_cells <- rbind(all_cells, data.frame(
        angle = angle,
        radius = xyl$r[1],
        x = x,
        y = y,
        id_layer = 20,
        id_cell = 1,
        order = 1.5,
        type = "xylem",
        id_group = i))
      i <-i+1
    }

    # protoxylem vessels are built on the outer stele rim
    protoxyl <- data.frame(r = max(all_cells$radius[all_cells$type == "stele"]) - (params$value[params$type == "cell_diameter" & params$name == "stele"])/2,
                           d = params$value[params$type == "cell_diameter" & params$name == "stele"])
    angle_seq_proto <- seq(from = 0, to = (2*pi), by = (2 * pi) / n_proto_xylem)
    for(angle in angle_seq_proto){
      x1 <- center + (protoxyl$r[1] * cos(angle))
      y1 <- center + (protoxyl$r[1] * sin(angle))
      #Find the closest stele cell and assign it as a protoxylem vessel
      all_cells <- all_cells %>%
        mutate(type = as.character(type)) %>%
        mutate(dist_protoxyl = sqrt((x-x1)^2 + (y-y1)^2)) %>%
        mutate(dist_protoxyl = ifelse(type == "stele", dist_protoxyl, 100)) %>%
        mutate(type = ifelse(dist_protoxyl == min(dist_protoxyl), "xylem", type))
    }

    # Phloem vessels are built between xylem ones
    phl <- data.frame(r = max(all_cells$radius[all_cells$type == "stele"]) - (params$value[params$type == "cell_diameter" & params$name == "stele"])/2,
                      d = params$value[params$type == "cell_diameter" & params$name == "stele"])
    angle_seq_ph <- seq(from = ((2 * pi) / n_proto_xylem ) /2, to = (2*pi), by = (2 * pi) / n_proto_xylem)
    for(angle in angle_seq_ph){
      x1 <- center + (phl$r[1] * cos(angle))
      y1 <- center + (phl$r[1] * sin(angle))
      #Find the closest stele cell and assign it as a phloem vessel
      all_cells <- all_cells %>%
        mutate(type = as.character(type)) %>%
        mutate(dist_phl = sqrt((x-x1)^2 + (y-y1)^2)) %>%
        mutate(dist_phl = ifelse(type == "stele", dist_phl, 100)) %>%
        mutate(type = ifelse(dist_phl == min(dist_phl), "phloem", type))

      # Get the compânion cells
      all_cells <- all_cells %>%
        mutate(type = as.character(type)) %>%
        mutate(dist_phl = sqrt((x-x1)^2 + (y-y1)^2)) %>%
        mutate(dist_phl = ifelse(type == "stele", dist_phl, 100)) %>%
        mutate(type = ifelse(dist_phl == min(dist_phl), "companion_cell", type))

      all_cells <- all_cells %>%
        mutate(type = as.character(type)) %>%
        mutate(dist_phl = sqrt((x-x1)^2 + (y-y1)^2)) %>%
        mutate(dist_phl = ifelse(type == "stele", dist_phl, 100)) %>%
        mutate(type = ifelse(dist_phl == min(dist_phl), "companion_cell", type))
    }

  }

  # Change the identity of stele cells to be replaced by xylem cells
  for(i in c(1:nrow(all_xylem))){
    # print(i)
    if(plant_type == 1){
      all_cells <- all_cells %>%
        mutate(type = ifelse(((x-all_xylem$x[i])^2 + (y - all_xylem$y[i])^2 < (all_xylem$d[i]/2)^2 & type == "stele"),
                             "xylem", type)) %>%
        mutate(id_group = ifelse(((x-all_xylem$x[i])^2 + (y - all_xylem$y[i])^2 < (all_xylem$d[i]/2)^2 & type == "xylem"),
                                 all_xylem$id_group[i], id_group))
    }else if(plant_type == 2){
      all_cells <- all_cells %>%
        filter(!((x-all_xylem$x[i])^2 + (y - all_xylem$y[i])^2 < (all_xylem$d[i]/1.5)^2 & type == "stele")) # find the cells inside the xylem poles and remove them
    }
  }

  if(plant_type == 2){
    all_xylem <- all_cells[all_cells$type == "xylem",]%>%
      mutate(radius = round(radius,2),
             d = c(xyl$d[1], rep(xyl$d[2:length(xyl$d)],
                                 max(unique(all_cells$id_group[all_cells$type == "xylem"]))-1)))%>%
      # Remove xylem cell that are to close to each other
      distinct(angle, radius, .keep_all = T)

    all_cells <- all_cells[all_cells$type != "xylem",]

    # Make circular frontier for xylem
    x_cir <- seq(-0.95,0.95,0.95/4)
    y_p <- sqrt(1-x_cir^2)
    y_m <- -sqrt(1-x_cir^2)
    xyl_frontier <- NULL
    k <- 1
    for (i_xyl in 1:nrow(all_xylem)) {
      tmp <- all_xylem[i_xyl,]
      cir <- tibble(x = rep(x_cir,2), y = c(y_p,y_m))
      cir <- cir*abs(tmp$d)/2
      cir$x <- cir$x + tmp$x
      cir$y <- cir$y + tmp$y
      xyl_frontier <- rbind(xyl_frontier, data.frame(angle = tmp$angle,
                                                     radius = tmp$radius,
                                                     x = cir$x,
                                                     y = cir$y,
                                                     id_layer = 20,
                                                     id_cell = 1,
                                                     type = "xylem",
                                                     order = 1.5,
                                                     id_group = k))
      k <- k + 1
    }
    all_cells <- rbind(all_cells, xyl_frontier)
  }

  angle_seq_passar <- round(seq(from = 0, to = (2*pi), by = (2 * pi) / n_proto_xylem),1)
  id_endo <- all_cells$id_cell[all_cells$type == "endodermis" & round(all_cells$angle,1) %in% angle_seq_passar]
  id_passar <- id_endo[sample(1:length(id_endo), n_passar, replace=T)]
  k <- 1
  while(n_passar != length(id_passar)){
    id_passar <- id_endo[sample(1:length(id_endo), n_passar, replace=T)]
    k <- k + 1
    if(k > 10){id_passar <- 1:n_passar}
  }
  all_cells$type[all_cells$id_cell %in% id_passar] <- "passage_cell"

  # reset the cell ids
  all_cells$id_cell <- c(1:nrow(all_cells))
  id_passar <- all_cells$id_cell[all_cells$type == "passage_cell"]

  t5 <- proc.time()




  #//////////////////////////////////////////////////////////////////////////////////////////////////
  # CREATE GEOMETRY ------
  if(verbatim) message("Creating the geometry")

  # Get the voronio data
  vtess <- deldir(all_cells$x, all_cells$y, digits = 8)
  if(is.null(vtess)){
    return(NULL)
  }

  # Remove the ouside cells, to get the voronoi data straight
  all_cells <- all_cells  %>%
    filter(type != "outside")
  test <- all_cells[all_cells$type == "xylem", ]


  # Get the size of the cells
  cell_size <- vtess$summary
  ids <- all_cells$id_cell
  cell_size$id_cell <- c(1:nrow(cell_size))
  all_cells <- merge(all_cells, cell_size[,c("id_cell", "dir.area")], by="id_cell")
  all_cells$area <- all_cells$dir.area
  all_cells$dist <- sqrt((all_cells$x - center)^2 + (all_cells$y - center)^2 )
  ids <- all_cells$id_cell

  sum(all_cells$area[all_cells$type == "xylem"])

  rs <- vtess$dirsgs[vtess$dirsgs$ind1 %in% ids |
                       vtess$dirsgs$ind2 %in% ids,]

  # Get the cooridnates for every nodes in the voronoi
  rs <- rs %>% arrange(ind1)
  rs2 <- data.frame(x = rs$x1, y=rs$y1, id_cell = rs$ind1)
  rs2 <- rbind(rs2, data.frame(x = rs$x2, y=rs$y2, id_cell = rs$ind1))
  rs2 <- rbind(rs2, data.frame(x = rs$x2, y=rs$y2, id_cell = rs$ind2))
  rs2 <- rbind(rs2, data.frame(x = rs$x1, y=rs$y1, id_cell = rs$ind2))

  rs2 <- merge(rs2, all_cells[,c("id_cell", "type", "area", "dist", "angle", "radius", "id_layer", "id_group")], by="id_cell")


  t6 <- proc.time()

  #//////////////////////////////////////////////////////////////////////////////////////////////////
  # CREATE AERENCHYMA -----
  if(verbatim) message("Killing cells to make aerenchyma")

  angle_inc <- (2 * pi) / n_aerenchyma_files
  angle_range_inc <- (2 * pi * proportion_aerenchyma / 100) / n_aerenchyma_files
  angle_range_inc_ini <- angle_range_inc
  safe_cortex_layer <- c(min(rs2$id_layer[rs2$type == "cortex"]))
  ini_cortex_area <- sum(all_cells$area[all_cells$type == "cortex" |
                                          all_cells$type == "endodermis" |
                                          all_cells$type == "exodermis" |
                                          all_cells$type == "epidermis"])
  surface_to_kill <- ini_cortex_area*proportion_aerenchyma
  cortex_area <- ini_cortex_area
  m_cortex_a <- mean(all_cells$area[all_cells$type == "cortex"])
  to_kill <- round(surface_to_kill/m_cortex_a)
  rs2$type <- as.character(rs2$type)
  cortex_layer <- sort(unique(rs2$id_layer[rs2$type == "cortex"]))
  `%!in%` <- compose(`!`, `%in%`)
  again <- 0
  try_nbr <- 0
  missing <- 0
  while(missing < ini_cortex_area - surface_to_kill + (ini_cortex_area - surface_to_kill)*0.1| again == 0){
    try_nbr <- try_nbr + 1
    angle <- runif(1, 0.6, 1) * pi/n_aerenchyma_files
    angle_range <- c(angle - angle_range_inc, angle + angle_range_inc)
    rs3 <- rs2
    for(j in c(1:n_aerenchyma_files)){
      angle_range <- c(angle - angle_range_inc, angle + angle_range_inc)
      rs3 <- rs3 %>%
        filter(!(id_layer %!in% safe_cortex_layer & type == "cortex" & angle > angle_range[1] & angle < angle_range[2]))
      angle <- angle + angle_inc
    }
    missing <- length(which(unique(all_cells$id_cell) %!in% unique(rs3$id_cell)))
    cortex_area <- ini_cortex_area - missing*m_cortex_a

    again <- 0
    if(missing > to_kill){
      angle_range_inc <- angle_range_inc - angle_range_inc*0.01
    }
    else if (missing < to_kill) {
      angle_range_inc <- angle_range_inc + angle_range_inc*0.01
    }else break

    if(cortex_area < (ini_cortex_area - surface_to_kill)-(ini_cortex_area - surface_to_kill)*0.1){
      again <- 1
    }

    if(try_nbr > 1000){
      warning("fail to kill cortex cells")
      return(NULL)
    }
  }

  rs2 <- rs3

  # angle <- runif(1, 0.6, 1) * pi/n_aerenchyma_files
  # for(j in c(1:n_aerenchyma_files)){
  #   angle_range <- c(angle - angle_range_inc, angle + angle_range_inc)
  #   rs2 <- rs2 %>%
  #     filter(!(id_layer != safe_cortex_layer & type == "cortex" & angle > angle_range[1] & angle < angle_range[2]))
  #   angle <- angle + angle_inc
  # }

  t7 <- proc.time()

  #//////////////////////////////////////////////////////////////////////////////////////////////////
  ## TIDY DATA ------

  if(verbatim) message("Tidying data before export")

  rs1 <- rs2 %>%
    dplyr::group_by(id_cell) %>%
    dplyr::mutate(my = mean(y)) %>%
    dplyr::mutate(mx = mean(x)) %>%
    dplyr::mutate(atan = atan2(y-my, x - mx)) %>%
    dplyr::arrange(id_cell, atan) %>%
    filter(!duplicated(atan))

  if(verbatim) message("Merging xylem vessels")
  groups <- NULL
  lost_points <- NULL
  for(i in c(1:max(rs1$id_group))){
    temp <- rs1 %>%
      filter(id_group == i) %>%
      group_by(x,y) %>%
      filter(row_number() == 1)

    if(nrow(temp) > 2){
      hull <- with(temp, ahull(x, y, alpha=1)) # Create a convex hull aroud the exiting xylem cells TODO

      # Get the points that were removed with the convexhull.
      # This will be used to track them later and moved them to the edge
      rem <- c(1:nrow(temp))
      rem <- rem[! rem %in% hull$arcs[,7]]
      temp2 <- temp[rem,]

      temp <- temp[hull$arcs[,7],] %>%
        mutate(id_cell = min(id_cell)) %>%
        dplyr::group_by(id_cell) %>%
        dplyr::mutate(area = sum(area))%>%
        dplyr::mutate(my = mean(y)) %>%
        dplyr::mutate(mx = mean(x)) %>%
        dplyr::mutate(atan = atan2(y-my, x - mx)) %>%
        arrange(id_cell, atan) %>%
        filter(!duplicated(atan))

      groups <- rbind(groups, temp)
      lost_points <- rbind(lost_points, temp2)
    }
  }
  rs1 <- rs1 %>%
    filter(id_group == 0)
  rs1 <- rbind(rs1, groups)


  # REMOVE THE WRONG XYLEM POINTS IN THE SURROIUNDING CELLS
  bad_points <- NULL # > TO MOVE
  not_so_bad_points <- NULL # > TO REMOVE
  for(i in c(1:nrow(lost_points))){
    if(nrow(rs1[rs1$x == lost_points$x[i],]) > 1){
      bad_points <- rbind(bad_points, (rs1[rs1$x == lost_points$x[i],])) # Get the points that need to be moved
    }
    if(nrow(rs1[rs1$x == lost_points$x[i],]) == 1){
      not_so_bad_points <- rbind(not_so_bad_points, (rs1[rs1$x == lost_points$x[i],])) # Get the points that need to be removed
    }
  }

  # Remove the "not so bad points"
  rs1 <- rs1 %>%
    mutate(remove = ifelse(x %in% not_so_bad_points$x & y %in% not_so_bad_points$y, "yes", "no")) %>%
    filter(remove == "no")

  # move bad points further away from closest centroid

  rs1 <- rs1 %>%
    mutate(bad = ifelse(x %in% bad_points$x & y %in% bad_points$y, "yes", "no"))
  groups <- groups%>%
    mutate(bad = ifelse(x %in% bad_points$x & y %in% bad_points$y, "yes", "no"))

  xylem_area <- groups%>%
    distinct(id_group, area, my, mx)%>%
    group_by(id_group)%>%
    summarise(area = sum(area))



  rs1 <- rs1 %>%
    dplyr::group_by(id_cell) %>%
    dplyr::mutate(my = mean(y)) %>%
    dplyr::mutate(mx = mean(x)) %>%
    dplyr::mutate(atan = atan2(y-my, x - mx)) %>%
    dplyr::arrange(id_cell, atan) %>%
    filter(!duplicated(atan))

  # bad_points[bad_points$x == out_cells$x[1],]

  t8 <- proc.time()

  # Reset the ids of the cells to be continuous
  ids <- data.frame(id_cell = unique(rs1$id_cell))
  ids$new <- c(1:nrow(ids))
  rs1 <- merge(rs1, ids, by="id_cell")
  rs1$id_cell <- rs1$new


  all_cells <- merge(all_cells, ids, by="id_cell")
  all_cells$id_cell <- all_cells$new


  tt <- proc.time()
  # outputing the inputs
  output <- data.frame(io = "input", name = params$name, type = params$type, value = params$value)

  head(rs1)
  for (i in unique(rs1$id_cell)) {
    tmp <- rs1[rs1$id_cell == i,]
    pol <- Polygon(tmp[, c("x","y")])
    rs1$area[rs1$id_cell == i] <-  pol@area
  }

  one_cells <- rs1%>%
    distinct(id_cell, type, id_group, area, .keep_all = TRUE)
  sum(one_cells$area[one_cells$type == "xylem"])
  all_cells <- merge(all_cells, one_cells, by = "id_cell")


  # adding the outputs by cell layers
  out <- ddply(all_cells, .(type.y), summarise, n_cells=length(type.y),
               layer_area = sum(area.y),
               cell_area = mean(area.y)) %>%
    mutate(name = type.y) %>%
    dplyr::select(-type.y) %>%
    gather(key = "type", value = "value", n_cells, layer_area, cell_area) %>%
    mutate(io = "output")%>%
    dplyr::select(io, everything())
  output <- rbind(output, out)

  time <- as.numeric(Sys.time()-t_1)

  # finaly we add the outputs for the whole section
  output <-rbind(output, data.frame(io="output", name="all", type="n_cells", value = nrow(all_cells)))

  output <-rbind(output, data.frame(io="output", name="stelar", type="layer_area",
                                    value = sum(all_cells$area.y[all_cells$order < 4])))
  output <-rbind(output, data.frame(io="output", name="cortex_to_epidermis", type="layer_area",
                                    value = sum(all_cells$area.y[all_cells$order > 3])))

  output <-rbind(output, data.frame(io="output", name="all", type="layer_area", value = sum(all_cells$area.y)))
  output <-rbind(output, data.frame(io="output", name="aerenchyma", type="layer_area", value = (ini_cortex_area - cortex_area)))
  output <-rbind(output, data.frame(io="output", name="aerenchyma", type="proportion", value = (ini_cortex_area - cortex_area)/ini_cortex_area))
  output <-rbind(output, data.frame(io="output", name="simulation", type="time", value = time))

  print(Sys.time()-t_1)

  rs1$sorting <- c(1:nrow(rs1))

  nodes <- rs1 %>%
    group_by(id_cell) %>%
    dplyr::mutate(xx = c(x[-1],x[1])) %>%
    dplyr::mutate(yy = c(y[-1],y[1]))

  nodes <- nodes %>%
    ungroup() %>%
    mutate(vertical = ifelse(x == xx, "true", "false")) %>%

    mutate(x1 = ifelse(x > xx, x, xx)) %>%
    mutate(x2 = ifelse(x > xx, xx, x)) %>%
    mutate(y1 = ifelse(x > xx, y, yy)) %>%
    mutate(y2 = ifelse(x > xx, yy, y)) %>%
    # Bug fix when wall is perfectly vertical
    mutate(y1 = ifelse(x == xx,
                       ifelse(y > yy, yy, y), y1)) %>%
    mutate(y2 = ifelse(x == xx,
                       ifelse(y > yy, y, yy), y2)) %>%
    mutate(wall_length = sqrt((x2-x1)^2 + (y2-y1)^2)) %>%
    mutate(wall_length2 = sqrt((xx-x)^2 + (yy-y)^2))


  if(maturity_x){
    tmp_m <- mean(nodes$area[nodes$type == "cortex"])
    nodes$type[nodes$type == "xylem" & nodes$area > tmp_m ] <- "stele"
  }


  walls <- nodes[!duplicated(nodes[,c('x1', 'x2', 'y1', 'y2')]),] %>%
    dplyr::select(c(x1, x2, y1, y2))

  walls$id_wall <- c(1:nrow(walls))

  nodes <- merge(nodes, walls, by=c("x1", "x2", "y1", "y2"))
  nodes <- nodes %>%
    arrange(sorting)

  t9 <- proc.time()


  if(verbatim){
    message("---------------")
    message("Time analysis: ")
    message("-- Loading")
    print(t2-t1)
    message("-- Cell layers")
    print(t3-t2)
    message("-- Cells")
    print(t4-t3)
    message("-- Xylem")
    print(t5-t4)
    message("-- Geometry")
    print(t6-t5)
    message("-- Aerenchyma")
    print(t7-t6)
    message("-- Merging")
    print(t8-t7)
    message("-- Tidying")
    print(t9-t8)
    message("------ All")
    print(t9-t1)
  }


  return(list(nodes = nodes,
              walls = walls,
              cells=all_cells,
              output = output))
}
