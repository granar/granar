
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

# the bad function that mess all the code
concavety <- function(data){

  nodes <- data%>%
    mutate(id_point = paste0(round(x,3),";", round(y,3)))
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
