
#' Create root anatomy
#'
#' This function creates a 2D root cross section anatomy based on global parameters. It is the core of GRANAR.
#' @param path Path to the XML file containing the different parameters for the simulation. Not needed if 'parameter' is set.  Default = NULL.
#' @param parameters Table with the different parameters. Not needed if 'path' is set.  Default = NULL.
#' @param verbatim Display textual information aboutt he simulation. Default = NULL.
#' @keywords root anatomy
#' @export
#' @examples
#' create_anatomy("PATH_TO_XLM_FILE")
#'


create_anatomy <- function(path = NULL,  # PAth
                           parameters = NULL,
                           verbatim = F){

  # Return NULL is no parameters are specified
  if( is.null(path) & is.null(parameters)){
    warning("Please specify a parameter set for the simulation")
    return(NULL)
  }

  if(verbatim) message("Loading parameters")

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
      to_find <- c("planttype", "randomness", "xylem", "phloem", "stele", "endodermis", "exodermis", "epidermis", "aerenchyma", "pericycle", "cortex")
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

  t1 <- proc.time()

  # PARAMETERS -----
  plant_type <- params$value[params$name == "planttype"]
  random_fact <- params$value[params$name == "randomness"] / 10 * params$value[params$name == "stele" & params$type == "cell_diameter"]
  stele_diameter <- params$value[params$name == "stele" & params$type == "layer_diameter"]
  n_aerenchyma_files <- params$value[params$name == "aerenchyma" & params$type == "n_files"]
  proportion_aerenchyma <- params$value[params$name == "aerenchyma" & params$type == "proportion"]
  n_xylem_files <- params$value[params$name == "xylem" & params$type == "n_files"]
  proto_meta_ratio <- params$value[params$name == "xylem" & params$type == "ratio"]
  n_proto_xylem <- round(n_xylem_files*proto_meta_ratio)

  t2 <- proc.time()


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
  layers$n_layers[layers$name == "stele"] <- round((stele_diameter/2) / layers$cell_diameter[layers$name == "stele"])
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
  for(i in c(2:nrow(all_layers))){
    # Update radius
    all_layers$radius[i] <- all_layers$radius[i-1] +
      all_layers$cell_diameter[i-1] / 2 +
      all_layers$cell_diameter[i] / 2
    if(all_layers$name[i] == "outside"){
      all_layers$radius[i] <- all_layers$radius[i-1] +
        all_layers$cell_diameter[i-1] / 2 +
        all_layers$cell_diameter[i] / 2
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

  ggplot(all_cells, aes(x, y, colour=type)) +
    geom_point() +
    coord_fixed()


  #//////////////////////////////////////////////////////////////////////////////////////////////////
  # CREATE XYLEM VESSELS -----
  # Create the xylem files
  # Get the extremes

  if(verbatim) message("Creating xylem and phloem vessels")

  all_cells$id_group <- 0

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
    xyl <- xyl %>% arrange(r)
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

    #modification 05/01
    r= max(all_cells$radius[all_cells$type == "stele"]) - (params$value[params$type == "cell_diameter" & params$name == "stele"])*2 -
      (params$value[params$type == "max_size" & params$name == "xylem"])/2
    xyl <- data.frame(r = r,
                      d = params$value[params$type == "max_size" & params$name == "xylem"])
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
            id_group = i
          )
        )
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
      phl <- data.frame(r = r,
                        d = params$value[params$type == "cell_diameter" & params$name == "stele"])
      angle_seq_ph <- seq(from = ((2 * pi) / n_xylem_files ) /2, to = (2*pi), by = (2 * pi) / n_xylem_files)
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
    all_cells <- all_cells %>%
       mutate(type = ifelse(((x-all_xylem$x[i])^2 + (y - all_xylem$y[i])^2 < (all_xylem$d[i]/2)^2 & type == "stele"),
                            "xylem", type)) %>%
      mutate(id_group = ifelse(((x-all_xylem$x[i])^2 + (y - all_xylem$y[i])^2 < (all_xylem$d[i]/2)^2 & type == "xylem"),
                               all_xylem$id_group[i], id_group)) #%>%

    # all_cells <- all_cells %>%
    # filter(!((x-all_xylem$x[i])^2 + (y - all_xylem$y[i])^2 < (all_xylem$d[i]/1.5)^2 & type == "stele")) # find the cells inside the xylem poles and remove them
  }


  # reset the cell ids
  all_cells$id_cell <- c(1:nrow(all_cells))

  t5 <- proc.time()



  #//////////////////////////////////////////////////////////////////////////////////////////////////
  # CREATE GEOMETRY ------
  if(verbatim) message("Creating the geometry")

  # Get the voronio data
  vtess <- deldir(all_cells$x, all_cells$y)

  # Remove the ouside cells, to get the voronoi data straight
  all_cells <- all_cells  %>%
    filter(type != "outside")

  # Get the size of the cells
  cell_size <- vtess$summary
  ids <- all_cells$id_cell
  cell_size$id_cell <- c(1:nrow(cell_size))
  all_cells <- merge(all_cells, cell_size[,c("id_cell", "dir.area")], by=c("id_cell"))
  all_cells$area <- all_cells$dir.area
  all_cells$dist <- sqrt((all_cells$x - center)^2 + (all_cells$y - center)^2 )


  ids <- all_cells$id_cell

  rs <- vtess$dirsgs[vtess$dirsgs$ind1 %in% ids |
                       vtess$dirsgs$ind2 %in% ids,]

  # Get the cooridnates for every nodes in the voronoi
  rs <- rs %>% arrange(ind1)
  rs2 <- data.frame(x = rs$x1, y=rs$y1, id_cell = rs$ind1)
  rs2 <- rbind(rs2, data.frame(x = rs$x2, y=rs$y2, id_cell = rs$ind1))
  rs2 <- rbind(rs2, data.frame(x = rs$x2, y=rs$y2, id_cell = rs$ind2))
  rs2 <- rbind(rs2, data.frame(x = rs$x1, y=rs$y1, id_cell = rs$ind2))
  rs2 <- rs2 %>% arrange(id_cell)
  rs2 <- merge(rs2, all_cells[,c("id_cell", "type", "area", "dist", "angle", "radius", "id_layer", "id_group")], by="id_cell")

  t6 <- proc.time()





  #//////////////////////////////////////////////////////////////////////////////////////////////////
  # CREATE AERENCHYMA -----
  if(verbatim) message("Killing cells to make aerenchyma")

  angle_inc <- (2 * pi) / n_aerenchyma_files
  angle_range_inc <- (2 * pi * proportion_aerenchyma / 2) / n_aerenchyma_files
  safe_cortex_layer <- min(rs2$id_layer[rs2$type == "cortex"])

  rs2$type <- as.character(rs2$type)
  angle <- runif(1, 0.6, 1) * pi/n_aerenchyma_files
  for(j in c(1:n_aerenchyma_files)){
    angle_range <- c(angle - angle_range_inc, angle + angle_range_inc)
    rs2 <- rs2 %>%
      filter(!(id_layer != safe_cortex_layer & type == "cortex" & angle > angle_range[1] & angle < angle_range[2]))
    angle <- angle + angle_inc
  }

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


  ggplot(rs1, aes(x, y, group=id_cell, fill=type)) +
    geom_polygon(colour="white")+
    coord_fixed()


  # Merge the adgecent xylem cells
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
        mutate(area = sum(area)) %>%
        mutate(id_cell = min(id_cell)) %>%
        dplyr::group_by(id_cell) %>%
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

 rs1%>%
    filter(type == "xylem" | type == "stele")%>%
    ggplot(aes(x,y))+
    geom_polygon(aes(group = id_cell, fill = type), colour="white")
 groups%>%
   filter(type == "xylem" | type == "stele")%>%
   ggplot(aes(x,y))+
   geom_polygon(aes(group = id_cell, fill = type), colour="white")+
   geom_point()+
   geom_point(aes(mx,my))

 # rs1 <- rs1%>%
 #  ungroup()
 # how_many_bad_points <- c(1:nrow(bad_points))
 # all_xyl <- rs1 %>% filter(type == "xylem")
 # pb = txtProgressBar(min = 0, max = length(how_many_bad_points), initial = 0, style = 3)
 #  for(i in how_many_bad_points){
 #    setTxtProgressBar(pb,i)
 #    bdp <- bad_points[i,]%>%
 #      ungroup()%>%
 #      mutate(bad = "yes")
#
 #   xyl <- groups %>%
 #     filter(bad == "yes")
 #   xyl2 <- rs1 %>%
 #     filter(id_cell %in% xyl$id_cell)  %>%
 #     filter(x %in% all_xyl$x & y %in% all_xyl$y)
#
 #    closest<-groups%>%
 #      ungroup()%>%
 #      filter(bad == "no")%>%
 #      mutate(x_to_p = sqrt((x-bdp$x)^2+(y-bdp$y)^2))%>%
 #      filter(x_to_p == min(x_to_p))
#
 #    # Neighbourhood
 #    right <- groups%>%
 #      ungroup()%>%
 #      filter(id_group == closest$id_group & bad == "no")%>%
 #      arrange(atan)
 #    right <- rbind(right[nrow(right), ], right, right[1,])
 #    tmp <- which(right$x == closest$x)
 #    ri <- right[tmp[1]+1, ]
 #    le <- right[tmp[length(tmp)]-1, ]
#
 #    tmp <- rbind(closest[,-15], ri, le, bdp)
 #    if(tmp$angle[4] >= 0 & tmp$angle[1] <= 0){print("BEPPPPP")}
 #    if(tmp$angle[4] <= 0 & tmp$angle[1] >= 0){print("BEPPPPP")}
 #    if(tmp$angle[4] - tmp$angle[1] > 0){
 #      closest_2 <- le
 #    }else{closest_2 <- ri}
#
 #    res <- line.line.intersection(c(out_cells$x[1], out_cells$y[1]),
 #                                  c(out_cells$x[2], out_cells$y[2]),
 #                                  c(xyl2$x[1], xyl2$y[1]),
 #                                  c(xyl2$x[2], xyl2$y[2]))
 #    # line <- lm(c(bdp$y, closest$my)~c(bdp$x, closest$mx))
 #    # if(is.na(line$coefficients[2])){line$coefficients[2] <- 0}
 #    # line2 <- lm(c(closest$y, closest_2$y)~c(closest$x, closest_2$x ))
 #    # if(is.na(line2$coefficients[2])){line2$coefficients[2] <- 0}
 #    # new_x <- unname((line2$coefficients[1]-line$coefficients[1])/
 #    #                   (line$coefficients[2]-line2$coefficients[2]))
 #    # new_y <- unname(line$coefficients[2]*new_x + line$coefficients[1])
#
 #     # Move the points in the stele cells
 #     bdp$x <- new_x
 #     bdp$y <- new_y
 #     bdp <- bdp[,-14]%>%
 #       mutate(remove = "no", bad = "new")
 #     ## Add a point into the xylem vessel
 #     new_xyl <- closest[,-c(14,15)]%>%
 #       mutate(remove = "new", bad = "new")
 #     new_xyl$x[1] <- new_x
 #     new_xyl$y[1] <- new_y
 #     rs1 <- rbind(rs1, new_xyl, bdp)
 #  }
#

#  all_xyl <- rs1 %>% filter(type == "xylem")
#   how_many_bad_points <- c(1:nrow(bad_points))
#   pb = txtProgressBar(min = 0, max = length(how_many_bad_points), initial = 0, style = 3)
#   for(i in how_many_bad_points){
#     setTxtProgressBar(pb,i)
#     # Find the closest xylem coordinate
#
#     xyl <- rs1 %>%
#       filter(x == bad_points$x[i] & y == bad_points$y[i])
#     xyl2 <- rs1 %>%
#       filter(id_cell %in% xyl$id_cell)  %>%
#       filter(x %in% all_xyl$x & y %in% all_xyl$y)
#    if(nrow(xyl) > 1 & nrow(xyl2) > 1){
#        xyl2 <- rs1 %>%
#          filter(id_cell %in% xyl$id_cell)  %>%
#          filter(x %in% all_xyl$x & y %in% all_xyl$y)
#
#       new_xyl <- all_xyl[all_xyl$x == xyl2$x[1],]
#
#       out_cells <- rs1 %>%
#         filter(id_cell %in% xyl$id_cell)
#
#       out_cells <- out_cells[duplicated(out_cells$x), ]
#       res <- line.line.intersection(c(out_cells$x[1], out_cells$y[1]),
#                                     c(out_cells$x[2], out_cells$y[2]),
#                                     c(xyl2$x[1], xyl2$y[1]),
#                                     c(xyl2$x[2], xyl2$y[2]))
#       res <- data.frame(x=res[1], y=res[2])
#
#       # Move the points in the stele cells
#       rs1$x[rs1$x == bad_points$x[i] & rs1$y == bad_points$y[i]] <- res$x
#       rs1$y[rs1$x == res$x & rs1$y == bad_points$y[i]] <- res$y
#
#       ## Add a point into the xylem vessel
#       new_xyl$x <- res$x
#       new_xyl$y <- res$y
#       rs1 <- rbind(rs1, new_xyl)
#     }else if(nrow(xyl) == 1){
#       rs1 <- rs1 %>%
#         mutate(remove = ifelse(x %in% xyl$x & y %in% xyl$y, "yes", "no")) %>%
#         filter(remove == "no")
#     }
#   }
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

  rs1%>%
    filter(type == "stele" | type == "xylem")%>%
    ggplot(aes(x,y))+
    geom_polygon(aes(group = id_cell, fill = type), colour = "white")+
    geom_point(aes(shape = bad), alpha = 0.5)+
    coord_fixed()

  all_cells <- merge(all_cells, ids, by="id_cell")
  all_cells$id_cell <- all_cells$new

  tt <- proc.time()
  # outputing the inputs
  output <- data.frame(io = "input", name = params$name, type = params$type, value = params$value)


  # adding the outputs by cell layers
  out <- ddply(all_cells, .(type), summarise, n_cells=length(type),
                                                layer_area = sum(area),
                                                cell_area = mean(area)) %>%
    mutate(name = type) %>%
    dplyr::select(-type) %>%
    gather(key = "type", value = "value", n_cells, layer_area, cell_area) %>%
    mutate(io = "output")%>%
    dplyr::select(io, everything())
  output <- rbind(output, out)

  one_cells <- rs1%>%
    distinct(id_cell, type, id_group, area, .keep_all = TRUE)
  type_area <- one_cells%>%
    dplyr::group_by(type)%>%
    dplyr::summarise(layer_area = sum(area))


  output$value[output$io == "output" & output$name == "xylem" & output$type == "layer_area"] <- type_area$layer_area[type_area$type == "xylem"]
  output$value[output$io == "output" & output$name == "stele" & output$type == "layer_area"] <- type_area$layer_area[type_area$type == "stele"]
  # finaly we add the outputs for the whole section
  output <-rbind(output, data.frame(io="output", name="all", type="n_cells", value = nrow(all_cells)))
  output <-rbind(output, data.frame(io="output", name="all", type="layer_area", value = sum(one_cells$area)))

  print(proc.time() - tt)

  rs1$sorting <- c(1:nrow(rs1))

  nodes <- rs1 %>%
    group_by(id_cell) %>%
    dplyr::mutate(xx = c(x[-1],x[1])) %>%
    dplyr::mutate(yy = c(y[-1],y[1]))
  #
  nodes %>%
    filter(type %in% c("stele", "xylem", "phloem")) %>%
    ggplot(aes(x, y, col=type, group=id_cell)) +
    geom_point() +
    geom_line()+
    coord_fixed()

  # plot(nodes2$xx, nodes$xx)

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
