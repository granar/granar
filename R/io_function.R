
# io_function for create_anatomy_3

`%!in%` <- compose(`!`, `%in%`)

cell_layer <- function(params){
  layers <- params %>%
    filter(type %in% c("cell_diameter","n_layers","order")) %>%
    spread(type, value) %>%
    filter(!is.na(n_layers)) %>%
    arrange(order)
  stele_diameter <- params$value[params$name == "stele" & params$type == "layer_diameter"]

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

  all_layers <- layer_info(all_layers)

  return(list(all_layers = all_layers , layers = layers))
}

layer_info <- function(all_layers){

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

  return(all_layers)
}

create_cells <- function(all_layers, random_fact){

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
    }else if(substr(all_layers$name[i], 1,6) == "cortex"){
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

  return(all_cells)
}

vascular <- function(all_cells, params, layers, center){


  n_xylem_files <- params$value[params$name == "xylem" & params$type == "n_files"]
  proto_meta_ratio <- params$value[params$name == "xylem" & params$type == "ratio"]
  n_proto_xylem <- round(n_xylem_files*proto_meta_ratio)
  plant_type <- params$value[params$name == "planttype"]
  if(length(all_cells$id_group[all_cells$type == "cortex"])> 0){
    k_max_cortex <- max(all_cells$id_group[all_cells$type == "cortex"])
  }else{k_max_cortex  <- 0}

  # length(all_cells$id_cell[all_cells$type == "inter_cellular_space"])

  all_cells%>%
    filter(id_group != 0)%>%
    ggplot()+
    geom_point(aes(x,y, colour = factor(id_group)))+coord_fixed()+guides(colour = F)

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
    # Phloem vessels are built between xylem ones
    phl <- data.frame(r = max(all_cells$radius[all_cells$type == "stele"]) - (params$value[params$type == "cell_diameter" & params$name == "stele"])/2,
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

      # Get the compagnion cells
      for (compa in 1:2) {
        all_cells <- all_cells %>%
          mutate(type = as.character(type)) %>%
          mutate(dist_phl = sqrt((x-x1)^2 + (y-y1)^2)) %>%
          mutate(dist_phl = ifelse(type == "stele", dist_phl, 100)) %>%
          mutate(type = ifelse(dist_phl == min(dist_phl), "companion_cell", type))
      }

      for (compa in 1:6) {
        all_cells <- all_cells %>%
          mutate(type = as.character(type)) %>%
          mutate(dist_phl = sqrt((x-x1)^2 + (y-y1)^2)) %>%
          mutate(dist_phl = ifelse(type == "stele", dist_phl, 100)) %>%
          mutate(type = ifelse(dist_phl == min(dist_phl), "cambium", type))
      }
      all_cells <- all_cells%>% select(-dist_phl)

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

    # remove stele cell inside metaxylem vessels
    for(i in c(1:nrow(all_xylem))){
      all_cells <- all_cells %>%
        filter(!((x-all_xylem$x[i])^2 + (y - all_xylem$y[i])^2 < (all_xylem$d[i]/1.5)^2 & type == "stele"))
    }

    # Make circular frontier for xylem
    x_cir <- seq(-0.95,0.95,0.95/4)
    y_p <- sqrt(1-x_cir^2)
    y_m <- -sqrt(1-x_cir^2)
    xyl_frontier <- NULL
    k <- 1
    for (i_xyl in 1:nrow(all_xylem)) {
      tmp <- all_xylem[i_xyl,]
      cir <- tibble(x = rep(x_cir,2), y = c(y_p,y_m))
      cir <- cir*abs(tmp$d*0.8)/2
      cir$x <- cir$x + tmp$x
      cir$y <- cir$y + tmp$y
      xyl_frontier <- rbind(xyl_frontier, data.frame(angle = tmp$angle,
                                                     radius = xyl$r[1],
                                                     x = cir$x,
                                                     y = cir$y,
                                                     id_layer = 20,
                                                     id_cell = 1,
                                                     type = "xylem",
                                                     order = 1.5,
                                                     id_group = k))

      k <- k + 1
    }
    xyl_frontier%>%
      ggplot()+
      geom_point(aes(x,y, colour = factor(id_group)))+
      geom_point(aes(x,y, colour = factor(id_group)), data = all_xylem)+
      coord_fixed()
    # add xyl frontier
    all_cells <- all_cells[all_cells$type != "xylem",]
    all_cells <- rbind(all_cells, xyl_frontier)
    all_cells$id_group[all_cells$type == "xylem"] <- all_cells$id_group[all_cells$type == "xylem" & all_cells$id_group != 0] + k_max_cortex

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

      # Get the compC"nion cells
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
    all_cells <- all_cells%>%select(-dist_phl, -dist_protoxyl)

  }

  # Change the identity of stele cells to be replaced by xylem cells
  if(plant_type == 2){
    # all_xylem <- all_xylem%>%
    #   filter(d < mean(all_xylem$d))
    # all_cells <- all_cells%>%
    #   filter(type != "xylem" | radius > mean(all_cells$radius[all_cells$type == "xylem"]) )
    #
    # all_xylem%>%
    #   ggplot()+
    #   geom_point(aes(x,y, colour = d))+
    #   coord_fixed()

    for(i in c(1:nrow(all_xylem))){
      # print(i)
      all_cells <- all_cells %>%
        filter(!((x-all_xylem$x[i])^2 + (y - all_xylem$y[i])^2 < (all_xylem$d[i]/1.5)^2 & type == "stele")) # find the cells inside the xylem poles and remove them
    }
  }

  if(plant_type == 2){
    #
    # xyl <- xyl[xyl$r < mean(xyl$r),]


    all_xylem <- all_cells[all_cells$type == "xylem",]%>%
      mutate(radius = round(radius,2),
             d = c(xyl$d[1], rep(xyl$d[2:length(xyl$d)],
                                 max(unique(all_cells$id_group[all_cells$type == "xylem"]))-1)))# %>%
    # Remove xylem cell that are to close to each other
    # filter(!duplicated(angle), !duplicated(radius))

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
      cir <- cir*abs(tmp$d*0.8)/2
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
      xyl_frontier%>%
        ggplot()+
        geom_point(aes(x,y, colour = factor(id_group)))+
        coord_fixed()
      k <- k + 1
    }
    all_cells <- rbind(all_cells, xyl_frontier)
    all_cells$id_group[all_cells$type == "xylem"] <- all_cells$id_group[all_cells$type == "xylem" & all_cells$id_group != 0] + k_max_cortex


    all_cells <- make_pith(all_cells, params, center)

  }
  all_cells %>%
    filter(id_group != 0)%>%
    ggplot()+
    geom_point(aes(x,y, colour = id_group))+
    coord_fixed()

  # reset the cell ids
  all_cells$id_cell <- c(1:nrow(all_cells))

  return(all_cells)
}

cell_voro <- function(all_cells, vtess, center){

  # Get the size of the cells
  #cell_size <- vtess$summary
  ids <- all_cells$id_cell
  #cell_size$id_cell <- c(1:nrow(cell_size))
  #all_cells <- merge(all_cells, cell_size[,c("id_cell", "dir.area")], by="id_cell")
  all_cells$area <- NA #all_cells$dir.area
  all_cells$dist <- sqrt((all_cells$x - center)^2 + (all_cells$y - center)^2 )

  rs <- vtess$dirsgs[vtess$dirsgs$ind1 %in% ids |
                       vtess$dirsgs$ind2 %in% ids,]

  # Get the cooridnates for every nodes in the voronoi
  rs <- rs %>% arrange(ind1)
  rs2 <- data.frame(x = rs$x1, y=rs$y1, id_cell = rs$ind1)
  rs2 <- rbind(rs2, data.frame(x = rs$x2, y=rs$y2, id_cell = rs$ind1))
  rs2 <- rbind(rs2, data.frame(x = rs$x2, y=rs$y2, id_cell = rs$ind2))
  rs2 <- rbind(rs2, data.frame(x = rs$x1, y=rs$y1, id_cell = rs$ind2))

  rs2 <- merge(rs2, all_cells[,c("id_cell", "type", "area", "dist", "angle", "radius", "id_layer", "id_group")], by="id_cell")
  rs2 %>%
    filter(type %in% c("cortex", "xylem", "inter_cellular_space"))%>%
    ggplot()+
    geom_point(aes(x,y, colour = factor(id_group), shape = type))+
    coord_fixed()+guides(colour = F)
  rs2 <- rs2%>%filter(type != "outside")

  return(list(all_cells = all_cells, rs2 = rs2))

}

smoothy_cells <- function(rs1){

  for(i in c(1:max(rs1$id_group))){
    temp <- rs1[rs1$id_group == i,]
    rs1 <- rs1[rs1$id_group != i,]
    if(nrow(temp)> 0){
      temp <- concavety(temp)
      rs1 <- rbind(rs1,temp)
    }
  }

  return(rs1)
}

concavety <- function(data){

  nodes <- data
  id_list <- NULL
  rock_roll <- NULL
  while(length(id_list) != length(unique(data$id_cell))){
    nodes <- nodes%>%
      filter(id_cell %!in% id_list)%>%
      dplyr::group_by(id_cell)%>%
      filter(!duplicated(id_point))

    double <- nodes%>%
      dplyr::group_by(id_point)%>%
      dplyr::summarise(n = n())%>%
      filter(n > 1)


    first <- nodes[nodes$y == max(nodes$y),][1,]# select a point in the outer border
    if(first$id_point %in% double$id_point){
      first <- nodes[nodes$x == max(nodes$x),][1,]
    }
    if(first$id_point %in% double$id_point){
      first <- nodes[nodes$x == min(nodes$x),][1,]
    }
    if(first$id_point %in% double$id_point){
      first <- nodes[nodes$y == min(nodes$y),][1,]
    }
    if(first$id_point %in% double$id_point){
      stop("huston, we have a problem")
    }

    first_point <- first$id_point[1] # full incremental should start somewhere
    first_cell <- first$id_cell[1]

    swing <- first_cell
    it <- T
    rock <- first[1,]
    cinp <- nodes[nodes$id_cell == swing,]
    last_pos <- which(cinp$id_point == first_point)-1
    if(length(last_pos) == 1){
      END_point <- ifelse(cinp$id_point[1] == first_point, last(cinp$id_point),
                          cinp$id_point[last_pos])
    }else{message = "multiple end point in the concave shape object"}
    while(it){
      if(is.na(last(rock$id_point))){print("breaks due to NA values in polygons coords")
        stop()}
      # junction point
      if(last(rock$id_point) %in% double$id_point){ # is it a juction point
        # unless they have only one point as a connection between the cells
        battle <- unique(nodes$id_cell[nodes$id_point == last(rock$id_point)])
        battle <- nodes[nodes$id_cell %in% battle,]
        fr <- battle%>%
          dplyr::group_by(id_point)%>%
          dplyr::summarise(n = n())%>%
          filter(n > 1)
        if (nrow(fr) == 1){ # only one point for contact between cells
          # do nothing
        }else{
          # print(swing)
          # switch from the previous id_cell to the next
          swing <- unique(nodes$id_cell[nodes$id_point == last(rock$id_point) &
                                          nodes$id_cell != rock$id_cell[rock$id_point == last(rock$id_point)]])


          # print(swing)
          if(length(swing) > 1){ # when triple point or quadri point
            # next cell should be the next in a anti-clockwise order
            depa <- nodes[nodes$id_cell %in% swing,]
            last_y <- nodes$y[nodes$id_point == last(rock$id_point)][1]
            last_x <- nodes$x[nodes$id_point == last(rock$id_point)][1]

            prev_rock <- rock[rock$id_point != last(rock$id_point),]
            av_la <- prev_rock[prev_rock$id_point == last(prev_rock$id_point),][1,]
            av_la$dist_point <- sqrt((av_la$x-last_x)^2+(av_la$y-last_y)^2)
            crit_angle <- ifelse(av_la$y-last_y >= 0, acos((av_la$x - last_x)/av_la$dist_point),
                                 2*pi-acos((av_la$x - last_x)/av_la$dist_point))

            depa <- depa%>%
              mutate(dist_point = sqrt((mx-last_x)^2+(my-last_y)^2),
                     crit = ifelse(my-last_y >= 0, acos((mx - last_x)/dist_point),
                                   2*pi-acos((mx - last_x)/dist_point)))
            swing <- ifelse(crit_angle < min(depa$crit), depa$id_cell[depa$crit == min(depa$crit)], #
                            ifelse(crit_angle < max(depa$crit), depa$id_cell[depa$crit == max(depa$crit)],
                                   depa$id_cell[depa$crit == min(depa$crit)] )
            )
          }
          # print(swing) # only one value
          yep <- nodes[nodes$id_cell == swing & nodes$id_point == last(rock$id_point),][1,] # take the next point
          rock <- rbind(rock, yep)
        }
      }

      cinp <- nodes[nodes$id_cell == swing,] # in this cell selection
      roll <- rbind(cinp[cinp$atan >= last(rock$atan),],cinp[cinp$atan < last(rock$atan),]) # sort points
      # print(roll)
      if(roll$id_point[1] == END_point){
        # print("will finish convave hull shortly")
        it <- F
        break()
      }
      pn <- roll[roll$id_point %!in% rock$id_point,][1,] # take the next point
      #print(pn)
      if(is.na(pn$x)){
        print(roll)
        pl <- nodes%>%
          ggplot()+
          geom_polygon(aes(x,y, group = id_cell, fill = id_cell),colour = "white", alpha = 0.5)+
          geom_point(aes(x,y), size = 2, alpha = 0.5, data = nodes[nodes$id_point %in% rock$id_point,])+

          coord_fixed()

        print(pl)
      }
      rock <- rbind(rock, pn)
      tata <- pn$atan


    }

    rock$id_cell <- first_cell
    rock$atan <- seq(-pi,pi, 2*pi/nrow(rock))[1:nrow(rock)]
    tmp_id_list <- unique(nodes$id_cell[which(point.in.polygon(nodes$x, nodes$y, rock$x,rock$y) > 0)])
    id_list <- c(id_list, tmp_id_list)
    rock_roll <- rbind(rock_roll, rock)

    for (i in unique(rock_roll$id_cell)) {
      tmp_check <- rock_roll[rock_roll$id_cell == i,]
      tmp_check <- tmp_check[!is.na(tmp_check$x), ]
      pol <- Polygon(tmp_check[, c("x","y")])
      rock_roll$area[rock_roll$id_cell == i] <-  pol@area
    }
    if(0 %in% rock_roll$area){
      message(paste0("id ",unique(data$id_cell), "have null area, try to fix bug"))
    }
    #
    # pl <- rock_roll%>%
    #   ggplot(aes(x,y))+
    #   geom_polygon(aes(x,y, group = id_cell, fill = factor(id_cell)), fill = "grey",alpha = 0.3,colour = "white", data = data)+
    #   geom_polygon(aes(x,y, group = id_cell, fill = factor(id_cell)),alpha = 0.6, colour = "white")+
    #   coord_fixed()+
    #   guides(fill = F)
    # print(pl)

  }


  return(rock_roll%>%select(colnames(data)))


}
