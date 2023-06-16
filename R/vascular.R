#' @title Add the vascular element in the stele
#'
#' remove stele cells and place vascular elements
#' add round shaped boundaries for xylem elements
#' @param params The input dataframe
#' @param all_cells The cellular dataframe
#' @param layers the layer dataframe
#' @param center The cross-section center
#' @keywords root
#' @export
#' @examples
#' all_cells <- vascular(all_cells, params, layers, center)
#'

vascular <- function(all_cells, params, layers, center){


  n_xylem_files <- params$value[params$name == "xylem" & params$type == "n_files"]
  proto_meta_ratio <- params$value[params$name == "xylem" & params$type == "ratio"]
  n_proto_xylem <- round(n_xylem_files*proto_meta_ratio)
  plant_type <- params$value[params$name == "planttype"]
  if(length(all_cells$id_group[all_cells$type == "cortex"])> 0){
    k_max_cortex <- max(all_cells$id_group[all_cells$type == "cortex"])
  }else{k_max_cortex  <- 0}

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
      r= max(all_cells$radius[all_cells$type == "stele"]) - (params$value[params$type == "cell_diameter" & params$name == "stele"])*1.5 -
        (params$value[params$type == "max_size" & params$name == "xylem"])/2
      xyl <- data.frame(r = r,
                        d = params$value[params$type == "max_size" & params$name == "xylem"])
    }
    all_xylem <- NULL
    angle_seq <- seq(from = 0, to = (2*pi)-(2 * pi) / n_xylem_files, by = (2 * pi) / n_xylem_files)
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

    # add xyl frontier
    all_cells <- all_cells[all_cells$type != "xylem",]
    all_cells <- rbind(all_cells, xyl_frontier)
    all_cells$id_group[all_cells$type == "xylem"] <- all_cells$id_group[all_cells$type == "xylem" & all_cells$id_group != 0] + k_max_cortex

    # protoxylem vessels are built on the outer stele rim
    protoxyl <- data.frame(r = max(all_cells$radius[all_cells$type == "stele"]) - (params$value[params$type == "cell_diameter" & params$name == "stele"])/2,
                           d = params$value[params$type == "cell_diameter" & params$name == "stele"])
    angle_seq_proto <- seq(from = 0, to = (2*pi)-(2 * pi) / n_proto_xylem, by = (2 * pi) / n_proto_xylem)
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
    for(i in c(1:nrow(all_xylem))){
      # print(i)
      all_cells <- all_cells %>%
        filter(!((x-all_xylem$x[i])^2 + (y - all_xylem$y[i])^2 < (all_xylem$d[i]/1.5)^2 & type == "stele")) # find the cells inside the xylem poles and remove them
    }
  }

  if(plant_type == 2){
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

      k <- k + 1
    }
    all_cells <- rbind(all_cells, xyl_frontier)
    all_cells$id_group[all_cells$type == "xylem"] <- all_cells$id_group[all_cells$type == "xylem" & all_cells$id_group != 0] + k_max_cortex

    if(length(params$value[params$name == "pith" & params$type == "layer_diameter"]) > 0){
      all_cells <- make_pith(all_cells, params, center)
    }
  }
  # reset the cell ids
  all_cells$id_cell <- c(1:nrow(all_cells))

  return(all_cells)
}
