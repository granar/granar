
`%!in%` <- compose(`!`, `%in%`)

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

# intercellular space

rondy_cortex <- function(params, all_cells, center){
  random_fact <- params$value[params$name == "randomness"] / 10
  cor_d <- params$value[params$name == "cortex" & params$type == "cell_diameter"]
  
  # calibration parameter
  # to_adjust1 <- params$value[params$name == "coefficient" & params$type == "icp_size"]
  # to_adjust2 <- params$value[params$name == "coefficient" & params$type == "icp_ratio"]
  
  all_cortex <- all_cells[all_cells$type %in% c("cortex","endodermis", "exodermis"),]
  
  icp_size <- params$value[params$name == "inter_cellular_space" & params$type == "size"]
  if(length(icp_size)> 0){
    scaling <- 1-(icp_size/cor_d)/2+0.07
    if(scaling >= 0.99){
      scaling = 0.99
    }
  }else{scaling <- 0.95}
  

  if(length(all_cells$id_group[all_cells$type == "xylem"]) > 0){
  k_max_xylem <- max(all_cells$id_group[all_cells$type == "xylem"])
  }else{k_max_xylem <- 0}

  ctess <- deldir(all_cortex$x, all_cortex$y, digits = 8)
  idc <- unique(all_cortex$id_cell)
  idc <- 1:length(idc)
  rc <- ctess$dirsgs[ctess$dirsgs$ind1 %in% idc |
                       ctess$dirsgs$ind2 %in% idc,]
  rc <- rc%>% arrange(ind1)
  rc2 <- data.frame(x = rc$x1, y=rc$y1, id_cell = rc$ind1)
  rc2 <- rbind(rc2, data.frame(x = rc$x2, y=rc$y2, id_cell = rc$ind1))
  rc2 <- rbind(rc2, data.frame(x = rc$x2, y=rc$y2, id_cell = rc$ind2))
  rc2 <- rbind(rc2, data.frame(x = rc$x1, y=rc$y1, id_cell = rc$ind2))
  
  inner <- min(all_cortex$radius[all_cortex$type %in% c("cortex")]) # [all_cortex$type %in% c("endodermis", "cortex")]
  outer <- max(all_cortex$radius[all_cortex$type %in% c("cortex")]) # no intercellular space between exo cortex and cortex

  rc2 <- rc2%>%mutate(euc = sqrt((x-center)^2+(y-center)^2))%>%
    dplyr::group_by(id_cell)%>%
    dplyr::mutate(mx = mean(x),
                  my = mean(y),
                  atan = atan2(y-my, x - mx))%>%
    dplyr::arrange(id_cell, atan)
  
  all_cortex$id_cell <- 1:nrow(all_cortex)
  rc1 <- merge(rc2, all_cortex[,c("id_cell", "type", "radius", "id_layer")], by="id_cell")
  
  rcin <- rc2%>%
    filter(euc > inner,
           euc < outer)%>%
    mutate(ID = paste0(round(x,2),round(y,2)))%>%
    filter(!duplicated(ID))
  
  rc1%>%
    ggplot()+
    geom_polygon(aes(x,y, group = id_cell, fill = type), colour = "white")+
    #geom_point(aes(x,y,colour = factor(id_cell)), size = 2, alpha = 0.3, data = all_cortex)+
    # geom_point(aes(x,y, colour = factor(id_group)), data = cor_frontier)+
    coord_fixed()+
    guides(colour = F)
  
  all_inter <- data.frame(angle = ifelse(rcin$y-center >= 0, acos((rcin$x - center)/rcin$euc),
                                         2*pi-acos((rcin$x - center)/rcin$euc)),
                          radius = rcin$euc,
                          x = rcin$x, y = rcin$y,
                          id_layer = all_cortex$id_layer[1]+0.5,
                          id_cell = 1:nrow(rcin),
                          type = "inter_cellular_space",
                          order = params$value[params$name == "cortex" & params$type == "order"]+0.5,
                          id_group = 0 # if too close, they should be merge but not now 
  )
  if(length(params$value[params$name =="inter_cellular_space"]) > 0){
    coef_icp <- (10 * params$value[params$name =="inter_cellular_space" & params$type == "ratio"]) # here to modulate icp proportion coeficient
    if(coef_icp > 1){coef_icp = 1}
    inter_cellular_proportion <- coef_icp*nrow(all_inter)
  }else{
    inter_cellular_proportion <- 0.5*nrow(all_inter)
  }
  if(inter_cellular_proportion == 0){
    all_inter <- NULL
  }else{
  to_keep <- sample(1:nrow(all_inter), round(inter_cellular_proportion), replace=F)
  all_inter <- all_inter[all_inter$id_cell %in% to_keep,]
  all_inter$id_point <- paste0(round(all_inter$x,3),";",round(all_inter$y,3)) # 1 5m precision
  all_inter <- all_inter%>%
    filter(!duplicated(id_point))%>%
    select(-id_point)
  all_inter%>%
    ggplot()+
    geom_point(aes(x,y))+
    coord_fixed()
  }

  nodes <- vertex(rc1%>%
                    filter(type == "cortex"))
  nodes <- nodes %>%
    filter(wall_length > 0)%>%
    mutate(m = (y2-y1)/(x2-x1),
           k = y1-m*x1,
           r_dist = abs(k+m*mx-my)/sqrt(1+m^2)) # distance between a point (mx,my) to a segment defined by to point (x1,y1; x2,y2)
  
  nodes%>%
    ggplot()+
    geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2))+
    coord_fixed()
  
  cor <- nodes%>%
    filter(wall_length > 0)%>%
    dplyr::group_by(id_cell)%>%
    dplyr::mutate(radius = min(r_dist))%>%
    filter(!duplicated(id_cell))
  
  cor$radius[cor$id_layer %in% c(inner, outer)] <- cor$radius[cor$id_layer %in% c(inner, outer)]*0.45
  
  cor$id_group = 1:nrow(cor)
  
  circus <- seq(-0.95,0.95,0.95/4)
  cir <- data.frame(x_cir = rep(circus,2*nrow(cor)))%>%
    mutate(y_cir = rep(c(sqrt(1-circus^2),-sqrt(1-circus^2)),nrow(cor)),
           # mx = rep(cor$mx,2*length(circus)),
           # my = rep(cor$my,2*length(circus)),
           id_group = sort(rep(1:nrow(cor),2*length(circus))))
  
  cor_frontier <- merge(cir, cor[,c("id_group", "radius", "mx", "my", "id_layer")], by = "id_group")%>%
    transmute(radius = radius,
              x = x_cir*radius*scaling+mx,
              y = y_cir*radius*scaling+my,
              euc = sqrt((mx-center)^2+(my-center)^2),
              angle = ifelse(my-center > 0,acos((mx - center)/euc),
                             2*pi-acos((mx - center)/euc)) ,
              id_layer = id_layer,
              id_cell = 1,
              type = "cortex",
              order = params$value[params$name == "cortex" & params$type == "order"],
              id_group = id_group
              )%>%
    select(-euc)
  
  cor_frontier%>%
    ggplot()+
    geom_point(aes(x,y, colour = factor(id_group)))+
    coord_fixed()+
    guides(colour = F)

  all_cells <- rbind(all_cells[all_cells$type != "cortex",], cor_frontier)
  all_cells$id_group[all_cells$type == "cortex"] <- all_cells$id_group[all_cells$type == "cortex" & all_cells$id_group != 0] + k_max_xylem 
  all_cells <- rbind(all_cells, all_inter)

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

create_aerenchyma <- function(all_cells, params, rs2){
  
  # rs2 <- rs1
  
  n_aerenchyma_files <- params$value[params$name == "aerenchyma" & params$type == "n_files"]
  proportion_aerenchyma <- params$value[params$name == "aerenchyma" & params$type == "proportion"]
  
  angle_inc <- (2 * pi) / n_aerenchyma_files
  angle_range_inc <- (2 * pi * proportion_aerenchyma / 100) / n_aerenchyma_files
  angle_range_inc_ini <- angle_range_inc
  safe_cortex_layer <- c(min(rs2$id_layer[rs2$type == "cortex"]))
  ini_cortex_area <- sum(all_cells$area[all_cells$type %in% c("cortex" ,"endodermis", "exodermis" , 
                                                              "epidermis", "inter_cellular_space")])
  surface_to_kill <- ini_cortex_area*proportion_aerenchyma
  cortex_area <- ini_cortex_area
  m_cortex_a <- mean(all_cells$area[all_cells$type == "cortex"])
  to_kill <- round(surface_to_kill/m_cortex_a)
  # rs2$type <- as.character(rs2$type)
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
    chyma <- NULL
    for(j in c(1:n_aerenchyma_files)){
      angle_range <- c(angle - angle_range_inc, angle + angle_range_inc)
      
      ma_cell <- rs3 %>%
        filter((id_layer %!in% safe_cortex_layer & type == "cortex" & angle > angle_range[1] & angle < angle_range[2]))
      rs3 <- rs3 %>%
        filter(!(id_layer %!in% safe_cortex_layer & type == "cortex" & angle > angle_range[1] & angle < angle_range[2]))
      
      if(nrow(ma_cell) > 0){
        chyma <- rbind(chyma, ma_cell%>%mutate(zone = j))
      }
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
  
  
  chyma$inside <- "out"
  for (i in 1:nrow(chyma)) {
    chyma$inside[i] <- ifelse(length(unique(chyma$id_cell[chyma$x == chyma$x[i] & chyma$zone == chyma$zone [i]])) == 3,"in", "out")
  }
  
  
  # angle <- runif(1, 0.6, 1) * pi/n_aerenchyma_files
  # for(j in c(1:n_aerenchyma_files)){
  #   angle_range <- c(angle - angle_range_inc, angle + angle_range_inc)
  #   rs2 <- rs2 %>%
  #     filter(!(id_layer != safe_cortex_layer & type == "cortex" & angle > angle_range[1] & angle < angle_range[2]))
  #   angle <- angle + angle_inc
  # }
  
  t7 <- proc.time()
  
  
  rs2 <- rs3
  
  return(list(chyma = chyma, rs2 = rs2, cortex_area = cortex_area))
}

fuzze_inter <- function(rs1){
  # saved_rs <- rs1
  space <- rs1[rs1$type == "inter_cellular_space",]
  if(nrow(space)> 0){
  
  space$id_point <- paste0(round(space$x,3), ";",round(space$y,3))
  double <- space%>%
    dplyr::group_by(id_point)%>%
    dplyr::summarise(n = n())%>%
    filter(n > 1)
  to_correct <- unique(space$id_cell[space$id_point %in% double$id_point])
  
  #  icpp <- unique(space$id_cell[space$id_point %in% double$id_point])
  if(nrow(double) > 0){
  done <- NULL
  itm <- 0
  for (btw in to_correct) {
    comu1 <- space$id_point[space$id_cell == btw & space$id_point %in% double$id_point]
    nei <- unique(space$id_cell[space$id_cell != btw & space$id_point %in% comu1])
    bou <- c(btw, nei)
    ke <- length(bou)
    te <- 0
    while(ke > te){
      te <- length(bou)
      comu1 <- space$id_point[space$id_cell %in% bou & space$id_point %in% double$id_point]
      nei <- unique(space$id_cell[space$id_cell %!in% bou & space$id_point %in% comu1])
      bou <- unique(c(bou, nei))
      ke <- length(bou) 
    }
    itm <- bou
    
    if(itm[1] %in% done){next()}
    # print(bou)
    tmp_cell <- space[space$id_cell %in% itm, ]
    for (i in unique(tmp_cell$id_cell)) {
      tmp <- tmp_cell[tmp_cell$id_cell == i,]
      tmp <- tmp[!is.na(tmp$x), ]
      pol <- Polygon(tmp[, c("x","y")])
      tmp_cell$area[tmp_cell$id_cell == i] <-  pol@area
    }

    tmp_cell <- tmp_cell[tmp_cell$area > 0, ]
    if(nrow(tmp_cell) > 0){
    tmp_cell <- tmp_cell%>%
      dplyr::group_by(id_cell)%>%
      dplyr::filter(!duplicated(id_point))
    
    # only works for merging two cells with a single slice in between
    four <- which(tmp_cell$id_point %in% double$id_point)
    if(length(four) == 4 ){
      seg <- tmp_cell[four,]
      tmp_cell$my <- mean(seg$y)
      tmp_cell$mx <- mean(seg$x)
      
      tmp_cell$id_cell <- itm[1]
      tmp_cell$atan = atan2(tmp_cell$y-tmp_cell$my, tmp_cell$x - tmp_cell$mx)
      tmp_cell <- tmp_cell%>%
        dplyr::arrange(atan)
      rs1 <- rs1[rs1$id_cell %!in% itm,]
      rs1 <- rbind(rs1,tmp_cell%>%select(-id_point))
      done <- c(done, itm)
      }else{
      # saved_cells <- tmp_cell
      # if(saved_cells$id_cell %in% c(6464,6462)){
      #   stop()
      # }
      # print(four)
      rs1 <- rs1[rs1$id_cell %!in% itm,]
      tmp_cell <- concavety(tmp_cell)
      rs1 <- rbind(rs1,tmp_cell%>%select(-id_point))
    }}else{
        rs1 <- rs1[rs1$id_cell %!in% itm,]
      }
      done <- c(done, itm)
    }
  # saved_cells%>%
  #   filter(type == "inter_cellular_space")%>%
  #   ggplot()+
  #   geom_polygon(aes(x,y, group = id_cell, fill = factor(id_cell)), colour = "white", alpha = 0.2)+
  #   coord_fixed()
  }
  
  # rs1%>%
  #   filter(type == "inter_cellular_space")%>%
  #   ggplot()+
  #   geom_polygon(aes(x,y, group = id_cell, fill = factor(id_cell)), colour = "white", alpha = 0.8)+
  #   coord_fixed()
  rs1 <- rs1[!is.na(rs1$x),]
  }
  
  return(rs1)
}

make_pith <- function(all_cells, params, center){
  if(params$value[params$name == "pith"][1] > 0){
    pith_size <- params$value[params$name == "pith" & params$type == "layer_diameter"]/2
    pcell <- params$value[params$name == "pith" & params$type == "cell_diameter"]
  }else{pith_size <- 0}
  
  if(pith_size > 0){
    
    
    xylem <- all_cells%>%
      filter(type == "xylem")
    
    xylem <- xylem%>%
      dplyr::group_by(id_group)%>%
      dplyr::mutate(mx = mean(x),
                    my = mean(y),
                    euc = sqrt((mx-center)^2+(my - center)^2))
    inner <- unique(xylem$id_group[xylem$euc < pith_size])
    all_cells <- all_cells[all_cells$type != "xylem" | all_cells$id_group %!in% inner,]
    

    n_pith_lay <- round(1+(pith_size-pcell/2)/pcell)
    pith_layer <- data.frame(name="stele",
               n_layers=rep(n_pith_lay, n_pith_lay),
               cell_diameter=pcell,
               order = 0.5)
    
    pith_layer <- layer_info(pith_layer)
    new_cells <- create_cells(all_layers = pith_layer, random_fact = 0.001)
    new_cells%>%
      ggplot()+geom_point(aes(x,y))+
      coord_fixed()
    new_center <- mean(new_cells$x[new_cells$angle == 0], new_cells$y[new_cells$angle == 0]) 
    
    new_cells$x <- new_cells$x-new_center+center
    new_cells$y <- new_cells$y-new_center+center
    new_cells$id_group <- 0
    
    all_cells <- all_cells[sqrt((all_cells$x-center)^2+(all_cells$y-center)^2) > pith_size,]
    all_cells <- rbind(new_cells, all_cells)
    all_cells$id_cell <- 1:nrow(all_cells)
    
    all_cells%>%
      # filter(id_group %!in% inner)%>%
      ggplot()+
      geom_point(aes(x,y, colour = type))+
      # geom_point(aes(x,y), colour = "red", alpha = 0.2, data =     xylem%>%
      #              filter(id_group %in% inner))+
      # geom_point(aes(x,y), colour = "green", data = new_cells)+
      coord_fixed()
  }

  return(all_cells)
}

root_hair <- function(rs1, params, center){
  random_fact <- params$value[params$name == "randomness"] / 10 * params$value[params$name == "stele" & params$type == "cell_diameter"]
  n_hair <- params$value[params$name == "hair" & params$type == "n_files"]
  len_hair <- params$value[params$name == "hair" & params$type == "length"]
  hair_r <- params$value[params$name == "epidermis" & params$type == "cell_diameter"]/2
  if(length(n_hair) != 0){
    x0 <- center
    y0 <- center
    id_h <- sample(unique(rs1$id_cell[rs1$type == "epidermis"]), n_hair)
    hairy <- NULL
    for (h in id_h) {
      # Select randomly which cell has a root hair
      tmp_cell <- rs1%>%filter(id_cell == h)
      tmp_cell <- tmp_cell%>%
        mutate(euc = sqrt((x0-x)^2 + (y0-y)^2),
               out = ifelse(euc > mean(euc), "out", "in"), # select point in the outer part of the cell
               rank = rank(euc))
      if(length(tmp_cell$out[tmp_cell$out == "out"])>2){
        junc <- tmp_cell[tmp_cell$out == "out" & tmp_cell$euc != max(tmp_cell$euc),]
        # draw line in between the center of the cross section and the cell mass center
        x1 <- tmp_cell$mx[1]
        y1 <- tmp_cell$my[1]
        m <- (y1-y0)/(x1-x0)
        d <- y0-m*x0
        # Length of the root hair
        l_hair <- len_hair+ len_hair*runif(1,-1,1)*random_fact*100
        sig <- (l_hair^2)*(1+m^2)-(y1-m*x1-d)^2
        sigy <- ((l_hair-hair_r)^2)*(1+m^2)-(y1-m*x1-d)^2
        x2_p <- (x1+y1*m-d*m+sqrt(sig))/(1+m^2)
        x2_m <- (x1+y1*m-d*m-sqrt(sig))/(1+m^2)
        y2_p <- (d+x1*m+y1*m^2+sqrt(sig)*m)/(1+m^2)
        y2_m <- (d+x1*m+y1*m^2-sqrt(sig)*m)/(1+m^2)
        di <- rbind(tibble(euc = sqrt((x0-x2_p)^2+(y0-y2_p)^2), pm = "p"), 
                    tibble(euc = sqrt((x0-x2_m)^2+(y0-y2_m)^2), pm = "m"))
        if(di$pm[di$euc == max(di$euc)] == "p"){
          tmp_cell$x[tmp_cell$euc == max(tmp_cell$euc)] <- x2_p
          tmp_cell$y[tmp_cell$euc == max(tmp_cell$euc)] <- y2_p
          
          #central point of the root hair apex
          x3 <- (x1+y1*m-d*m+sqrt(sigy))/(1+m^2)
          y3 <- (d+x1*m+y1*m^2+sqrt(sigy)*m)/(1+m^2)
          # perpendicular line
          inv_m <- -1/m
          d3 <- y3-inv_m*x3
          # distance
          s <- (hair_r^2)*(1+inv_m^2)-(y3-inv_m*x3-d3)^2
          #coordinate of the point
          x4 <- (x3+y3*inv_m-d3*inv_m+sqrt(s))/(1+inv_m^2)
          x5 <- (x3+y3*inv_m-d3*inv_m-sqrt(s))/(1+inv_m^2)
          y4 <- (d3+x3*inv_m+y3*inv_m^2+sqrt(s)*inv_m)/(1+inv_m^2)
          y5 <- (d3+x3*inv_m+y3*inv_m^2-sqrt(s)*inv_m)/(1+inv_m^2)
          
          junc$x[junc$atan == min(junc$atan)] <- x4
          junc$y[junc$atan == min(junc$atan)] <- y4
          junc$x[junc$atan == max(junc$atan)] <- x5
          junc$y[junc$atan == max(junc$atan)] <- y5
          
          tmp_cell <- rbind(tmp_cell, junc)%>%
            mutate(my = mean(y),
                   mx = mean(x),
                   atan = atan2(y-my, x - mx))%>%
            arrange(atan)%>%
            filter(!duplicated(atan))
        }else{
          tmp_cell$x[tmp_cell$euc == max(tmp_cell$euc)] <- x2_m
          tmp_cell$y[tmp_cell$euc == max(tmp_cell$euc)] <- y2_m
          
          #central point of the root hair apex
          x3 <- (x1+y1*m-d*m-sqrt(sigy))/(1+m^2)
          y3 <- (d+x1*m+y1*m^2-sqrt(sigy)*m)/(1+m^2)
          # perpendicular line
          inv_m <- -1/m
          d3 <- y3-inv_m*x3
          # distance
          s <- (hair_r^2)*(1+inv_m^2)-(y3-inv_m*x3-d3)^2
          #coordinate of the point
          x4 <- (x3+y3*inv_m-d3*inv_m+sqrt(s))/(1+inv_m^2)
          x5 <- (x3+y3*inv_m-d3*inv_m-sqrt(s))/(1+inv_m^2)
          y4 <- (d3+x3*inv_m+y3*inv_m^2+sqrt(s)*inv_m)/(1+inv_m^2)
          y5 <- (d3+x3*inv_m+y3*inv_m^2-sqrt(s)*inv_m)/(1+inv_m^2)
          
          junc$x[junc$atan == min(junc$atan)] <- x4
          junc$y[junc$atan == min(junc$atan)] <- y4
          junc$x[junc$atan == max(junc$atan)] <- x5
          junc$y[junc$atan == max(junc$atan)] <- y5
          
          tmp_cell <- rbind(tmp_cell, junc)%>%
            mutate(my = mean(y),
                   mx = mean(x),
                   atan = atan2(y-my, x - mx))%>%
            arrange(atan)%>%
            filter(!duplicated(atan))
          
        }
        
      }
      # make a binding table with all root hair
      hairy <- rbind(hairy, tmp_cell)

    }
    rs1 <- rs1[rs1$id_cell %!in% id_h,]
    rs1$euc <- sqrt((x0-rs1$x)^2 + (y0-rs1$y)^2)
    rs1$out <- "none"
    rs8 <- rbind(rs1,hairy%>%select(-rank))%>%
      arrange(id_cell, atan)
  }
  
  return(rs8)
}

smoothy_cells <- function(rs1){
  # saved_rs <- rs1
  # rs1 <- saved_rs
  rs1$id_point <- paste0(round(rs1$x,3),";",round(rs1$y,3))
  rs1%>%
    ggplot(aes(x, y, group=id_cell, fill = factor(id_group)))+
    geom_polygon(colour="white")+
    coord_fixed()+
    guides(fill = F)
  
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
  

vertex <- function(rs1){
  nodes <- rs1 %>%
    mutate(id_point= paste0(round(rs1$x,3),";",round(rs1$y,3)))
    group_by(id_cell) %>%
    filter(!duplicatded(id_point))%>%
    dplyr::mutate(xx = c(x[-1],x[1])) %>%
    dplyr::mutate(yy = c(y[-1],y[1]))%>%
      select(-id_point)
  
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
    mutate(wall_length2 = sqrt((xx-x)^2 + (yy-y)^2),
           slope = (y2-y1)/(x2-x1),
           intercept = y1 - slope*x1)
  return(nodes)
}

concavety <- function(data){

  # data %>%
  #   ggplot()+
  #   geom_polygon(aes(x,y, group = id_cell, fill = id_cell), alpha = 0.2,colour = "white")+
  #   geom_point(aes(x,y), data = rock_roll)+
  #   #geom_polygon(aes(x,y),fill = "blue",alpha = 0.5, colour = "white", data = rock_roll)+
  #   coord_fixed()
  #  data <- saved_cells
  
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

aerenchyma <- function(params, rs1){
  
  # saved_rs <- rs1
  # rs1 <- saved_rs
  
  aer_type <- params$value[params$name == "aerenchyma" & params$type == "type"]
  if(length(aer_type)== 0){
    aer_type <- params$value[params$name == "planttype" & params$type == "param"]
  }
  n_aerenchyma_files <- params$value[params$name == "aerenchyma" & params$type == "n_files"]
  proportion_aerenchyma <- params$value[params$name == "aerenchyma" & params$type == "proportion"]
  
  angle_inc <- (2 * pi) / n_aerenchyma_files
  safe_cortex_layer <- c(min(rs1$id_layer[rs1$type == "cortex"]))
  last_cortex_layer <- c(max(rs1$id_layer[rs1$type == "cortex"]))
  area_all <- rs1%>%
    filter(type %in% c("cortex", "inter_cellular_space", "exodermis", "epidermis"))%>%
    dplyr::group_by(id_cell)%>%
    dplyr::summarise(area = mean(area))
  ini_cortex_area <- sum(area_all$area)
  surface_to_kill <- ini_cortex_area*proportion_aerenchyma
  stk_zone = surface_to_kill/n_aerenchyma_files
  # (2 * pi * proportion_aerenchyma / 30) / n_aerenchyma_files
  if (aer_type == 1){
  small_r <- mean(rs1$dist[rs1$id_layer == safe_cortex_layer])+0.5*mean(rs1$radius[rs1$id_layer == safe_cortex_layer])
  big_R <- mean(rs1$dist[rs1$id_layer == last_cortex_layer])+0.5*mean(rs1$radius[rs1$id_layer == last_cortex_layer])
  angle_range_inc <- stk_zone/(big_R^2-small_r^2)
  
  }
  if(aer_type == 2){
    angle_range_inc <- (2 * pi * proportion_aerenchyma / 100) / n_aerenchyma_files
  }
  
  angle_range_inc_ini <- angle_range_inc
  cortex_area <- ini_cortex_area
  angle <- runif(1, 0.6, 1) * pi/n_aerenchyma_files
  angle_range <- c(angle - angle_range_inc, angle + angle_range_inc)
  
  rs1$id_point <- paste0(round(rs1$x,3), ";", round(rs1$y,3))
  id_group_max <- max(rs1$id_group)
  
  
  for(j in c(1:n_aerenchyma_files)){
    # saved_rs1 <- rs1
    possi <- rs1[rs1$id_layer %!in% safe_cortex_layer & rs1$type %in% c("cortex", "inter_cellular_space"),]
    gotogo <- T
    n_try <- 0
    while(gotogo){
    angle_range <- c(angle - angle_range_inc, angle + angle_range_inc)
    ma_cell <- rs1 %>%
      filter(id_layer %!in% safe_cortex_layer ,
             type == "cortex" ,
             angle > angle_range[1] ,
             angle < angle_range[2])
    are <- sum_area(ma_cell)
      if(are > stk_zone){
        angle_range_inc <- angle_range_inc*0.9
        n_try <- n_try + 1
        if(n_try > 100){
          gotogo <- F
        }
        next()
      }
      if(nrow(ma_cell) > 0){
      while(are <= stk_zone ){
        are <- sum_area(ma_cell)
        pointy <- unique(ma_cell$id_point)
        done <- unique(ma_cell$id_cell) # id from potential cell to be named aer
        nei <- unique(possi$id_cell[possi$id_cell %!in% done & possi$id_point %in% pointy]) # neigh cells and inter_cell
        icp <- unique(possi$id_cell[possi$id_cell %in% nei & possi$type == "inter_cellular_space"]) # neigh inter_cell
        pre <- unique(possi$id_cell[possi$id_cell %in% c(done,icp)]) # merge potential and inter_cell
        
        done <- c(done, nei) # full potential
        ma_cell <- possi[possi$id_cell %in% done,]
        are <- sum_area(ma_cell)
        
        # pl <- ggplot()+
        #   geom_polygon(aes(x,y,group = id_cell, fill = id_cell), colour = "white", alpha = 0.8,data = ma_cell)+
        #   geom_polygon(aes(x,y,group = id_cell, fill = id_cell), colour = "white", alpha = 0.1,data = rs1)+
        #   geom_polygon(aes(x,y,group = id_cell), fill = "red",colour = "white", alpha = 0.1,data = rs1%>%
        #                  filter(id_layer == safe_cortex_layer))+
        #   coord_fixed()
        # print(pl)
        
        n_try <- n_try + 1
        if(n_try > 100){
          are = stk_zone
        }
      }
        # once it overflow the threshold limit
        ma_cell <- possi[possi$id_cell %in% pre,] # this default potential cell and inter_cell are selected
        are_b <- sum_area(ma_cell)
        # adjust to include the right amount of cells as aerenchyma
        if(are_b > stk_zone){to_be_added <- NULL}else{
        
        next_area <- possi[possi$id_cell %in% done & possi$id_cell %!in% pre,]
        alm <- next_area%>%
          dplyr::group_by(id_cell)%>%
          dplyr::summarise(area = mean(area))%>%
          mutate(sarea = cumsum(area),
                 near = abs(sarea+are_b-stk_zone))
        jk <- which(alm$near == min(alm$near))
        
        to_be_added <- alm$id_cell[1:jk]
        ma_cell <- possi[possi$id_cell %in% c(pre,to_be_added),]
        pointy <- unique(ma_cell$id_point)
        done <- unique(ma_cell$id_cell) # id from potential cell to be named aer
        nei <- unique(possi$id_cell[possi$id_cell %!in% done & possi$id_point %in% pointy]) # neigh cells and inter_cell
        icp <- unique(possi$id_cell[possi$id_cell %in% nei & possi$type == "inter_cellular_space"])
        pre <- unique(possi$id_cell[possi$id_cell %in% c(done,icp)])
        ma_cell <- possi[possi$id_cell %in% pre,]
        }
        
        rs1 <- rs1[rs1$id_cell %!in% pre,]
        gotogo <- F # we are done for this part
        #ma_cell = saved_ma_cell <- ma_cell
        aer <- concavety(ma_cell)
        aer$type <- as.character("aerenchyma")
        aer$id_group <- j+id_group_max
        
        
        # print(unique(rs1$id_cell[rs1$id_cell %in% pre]))
        # print(str(aer))
        # print(str(rs1))
        rs1 <- rbind(rs1, aer)
        
        # print(aer%>%
        #   ggplot(aes(x,y))+
        #   geom_polygon(aes(x,y,group = id_cell, fill = type), colour = "white", alpha = 0.8,data = rs1)+
        #   geom_polygon(colour = "white", fill = "red")+
        #   coord_fixed())
        # 
      }else{
        angle_range_inc <- angle_range_inc*1.1
        n_try <- n_try + 1
        if(n_try > 100){
          gotogo <- F
        }
      }
      
    }
    angle <- angle + angle_inc
    message(paste0(j,"/",n_aerenchyma_files))
    # j <- j +1
  }
  rs <- rs1
  
  return(rs)
  
}

sum_area <- function(cells){
  area_tmp <- cells %>%
    dplyr::group_by(id_cell)%>%
    dplyr::summarise(area = mean(area))%>%
    arrange(area, decreasing = T)
  are <- sum(area_tmp$area)
  return(are)
}

septa <- function(rs1){
  data <- rs1
  data$id_point <- paste0(round(data$x,3),";",round(data$y,3))
  data$neib1 = data$neib2 <- NA
  data$neib1_type = data$neib2_type <- NA

  for(i in which(data$type == "aerenchyma")){
    ne <- unique(data$id_cell[data$id_point == data$id_point[i] & data$id_cell != data$id_cell[i]])
    data$neib1[i] <- ne[1]
    data$neib1_type[i] <- data$type[data$id_cell == ne[1]][1]
    data$neib2[i] <- ne[2]
    data$neib2_type[i] <- data$type[data$id_cell == ne[2]][1]
    if(!is.na(ne[3])){message("quadri point")
      print(ne)}
  }
  # all triple points and all points next to cortex cells
  must <- data[!is.na(data$neib2),]
  must <- rbind(must, data[data$neib1_type != "aerenchyma",])
  
  # add some noise 
  noise <- data[!(data$neib1_type != "aerenchyma" | !is.na(data$neib2)),]
  noise_point <- unique(noise$id_point)
  n_noise <- round(length(noise_point)*0.05)
  noise_keep <- sample(noise_point, n_noise, replace=F)
  noise <- noise[noise$id_point %in% noise_keep,]
  must <- rbind(must, noise)
  
  data <- data[data$type != "aerenchyma",]
  data <- rbind(data,must)%>%
    arrange(id_cell,atan)
  
  
  data %>%
    filter(type == "aerenchyma",
           !is.na(neib1))%>%
    ggplot()+
    geom_polygon(aes(x,y, group = id_cell, fill = type),alpha = 0.6, colour = "white", data = rs1)+
    geom_polygon(aes(x,y, group = id_cell, fill = type), colour = "white")+
    geom_point(aes(x,y), colour = "red", data = data[!is.na(data$neib2),])+
    geom_point(aes(x,y), colour = "blue", data = noise)+
    coord_fixed()
  
  
  return(data%>%select(colnames(rs1)))
}


pv_ready <- function(rs1){
  data <- rs1
  
  data$id_point <- paste0(round(data$x,3),";",round(data$y,3))
  data$neib1 = data$neib2 = data$neib3 <- NA
  data$neib1_type = data$neib2_type = data$neib3_type <- NA
  
  for(i in 1:nrow(data)){
    ne <- unique(data$id_cell[data$id_point == data$id_point[i] & data$id_cell != data$id_cell[i]])
    data$neib1[i] <- ne[1]
    data$neib1_type[i] <- data$type[data$id_cell == ne[1]][1]
    data$neib2[i] <- ne[2]
    data$neib2_type[i] <- data$type[data$id_cell == ne[2]][1]
    data$neib3[i] <- ne[3]
    data$neib3_type[i] <- data$type[data$id_cell == ne[3]][1]
  }
  
  junc_point <- data[!is.na(data$neib2) | (data$type == "epidermis" & !is.na(data$neib1)),]
  
  data %>%
    ggplot()+
    geom_polygon(aes(x,y,fill = type, group = id_cell), colour = "white")+
    geom_point(aes(x,y), colour = "red", data = junc_point)+
    geom_point(aes(x,y), colour = "blue", data = data[data$id_point %!in% junc_point$id_point,])+
    coord_fixed()+
    theme_classic()
  
  nodes <- data
  
  nodes<- nodes%>%
    dplyr::group_by(id_cell)%>%
    dplyr::filter(!duplicated(id_point))%>%
    dplyr::mutate(x1 = x,
                  y1 = y,
                  x2 = c(x[-1],x[1]),
                  y2 = c(y[-1],y[1]),
                  id_point2 = paste0(round(x2,3),";",round(y2,3)),
                  w_length = sqrt((x1-x2)^2+(y1-y2)^2)) # %>% #hard coded
    # dplyr::mutate(x3 = ifelse(id_point2 %in% junc_point$id_point | id_point2 == "NA;NA",NA,c(x2[-1],x2[1])),
    #               y3 = ifelse(id_point2 %in% junc_point$id_point | id_point2 == "NA;NA",NA,c(y2[-1],y2[1])),
    #               id_point3 = paste0(round(x3,3),";",round(y3,3)),
    #               x4 = ifelse(id_point3 %in% junc_point$id_point | id_point3 == "NA;NA",NA,c(x3[-1],x3[1])),
    #               y4 = ifelse(id_point3 %in% junc_point$id_point | id_point3 == "NA;NA",NA,c(y3[-1],y3[1])),
    #               id_point4 = paste0(round(x4,4),";",round(y4,4)),
    #               x5 = ifelse(id_point4 %in% junc_point$id_point | id_point4 == "NA;NA",NA,c(x4[-1],x4[1])),
    #               y5 = ifelse(id_point4 %in% junc_point$id_point | id_point4 == "NA;NA",NA,c(y4[-1],y4[1])),
    #               id_point5 = paste0(round(x5,5),";",round(y5,5)),
    #               x6 = ifelse(id_point5 %in% junc_point$id_point | id_point5 == "NA;NA",NA,c(x5[-1],x5[1])),
    #               y6 = ifelse(id_point5 %in% junc_point$id_point | id_point5 == "NA;NA",NA,c(y5[-1],y5[1])),
    #               id_point6 = paste0(round(x6,6),";",round(y6,6)),
    #               x7 = ifelse(id_point6 %in% junc_point$id_point | id_point6 == "NA;NA",NA,c(x6[-1],x6[1])),
    #               y7 = ifelse(id_point6 %in% junc_point$id_point | id_point6 == "NA;NA",NA,c(y6[-1],y6[1])),
    #               id_point7 = paste0(round(x7,7),";",round(y7,7)),
    #               x8 = ifelse(id_point7 %in% junc_point$id_point | id_point7 == "NA;NA",NA,c(x7[-1],x7[1])),
    #               y8 = ifelse(id_point7 %in% junc_point$id_point | id_point7 == "NA;NA",NA,c(y7[-1],y7[1])),
    #               id_point8 = paste0(round(x8,8),";",round(y8,8)),
    #               x9 = ifelse(id_point8 %in% junc_point$id_point | id_point8 == "NA;NA",NA,c(x8[-1],x8[1])),
    #               y9 = ifelse(id_point8 %in% junc_point$id_point | id_point8 == "NA;NA",NA,c(y8[-1],y8[1])),
    #               id_point9 = paste0(round(x9,9),";",round(y9,9))
    #               )%>%
    # dplyr::filter(id_point %in% junc_point$id_point,
    #               w_length > 0)# only one wall per junction point
  
  on_go <- T
  while(on_go){
    
    last_point <- nodes%>% 
      select(id_cell, last(starts_with("id_point")),
             last(starts_with("x")),
             last(starts_with("y"))) # take last colmn
    n <- colnames( last_point)
    n <- unique(parse_number(n[-1])) # get the incremental value
    colnames(last_point) <- c("id_cell","id_point", "x","y") # generic col names
    
    last_point <- last_point%>%
      dplyr::group_by(id_cell)%>%
      dplyr::mutate(xx = ifelse(id_point %in% junc_point$id_point | id_point == "NA;NA",NA,c(x[-1],x[1])),
             yy = ifelse(id_point %in% junc_point$id_point | id_point == "NA;NA",NA,c(y[-1],y[1])),
             id_pointxy = paste0(round(xx,3),";",round(yy,3)))%>%
      ungroup()%>%
      select(xx,yy,id_pointxy)
    
    on_go <- length(which(last_point$id_pointxy != "NA;NA"))>0
    c_name <- paste0(t(c("x", "y", "id_point")),n+1)
    colnames(last_point) <- c(c_name)

    nodes <- cbind(nodes%>%
                     ungroup(),last_point)
  }
  
  nodes <- nodes%>%
    dplyr::filter(id_point %in% junc_point$id_point,
                  w_length > 0)# only one wall per junction point
  
  nodes%>%
    ggplot()+
    geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2))+
    geom_segment(aes(x = x2, xend = x3, y = y2, yend = y3))+
    geom_segment(aes(x = x3, xend = x4, y = y3, yend = y4))+
    coord_fixed()
  
  
  nodus <- NULL
  more <- nodes 
  k <- 2
  while(nrow(more) > 0){
    print(k)
    last_x <- paste0("x",k) # the last point has coord_x 
    last_y <- paste0("y",k) # the last point has coord_y 
    tag_x <- paste0("x",k+1)  # the next possible points
    
    tmp <- more[is.na(more[,tag_x]),] # if next point is "nan" then  
    tmp <- tmp%>%
      # dplyr::group_by(id_cell) %>%
      dplyr::mutate(xx = tmp[,last_x])%>%
      dplyr::mutate(yy = tmp[,last_y])
    
    deto <- tmp%>%select(-x1,-y1,-c(last_x),-c(last_y)) # remove points
    h <- 1
    # dealing with first and last point
      tmp_h <- tmp%>%
        mutate(xh = ifelse(x > xx, x, xx)) %>%
        mutate(yh = ifelse(x > xx, y, yy)) %>%
        mutate(yh = ifelse(x == xx, ifelse (y > yy, yy, y), yh))%>%
        mutate(xlh = ifelse(x > xx, xx, x)) %>%
        mutate(ylh = ifelse(x > xx, yy, y)) %>%
        mutate(ylh = ifelse(x == xx, ifelse (y > yy, y, yy), ylh))%>%
        select(xh,yh, xlh, ylh)
      colnames(tmp_h) <- c(paste0("x",h),paste0("y",h), paste0("x",k),paste0("y",k))
      print(c(paste0("x",h),paste0("y",h), paste0("x",k),paste0("y",k)))
      deto <- cbind(deto, tmp_h)
    if(k >= 4){
      for(h in c(2:floor(k/2))){
        h_x <- paste0("x",h)
        lh_x <- paste0("x",k+1-h)
        h_y <- paste0("y",h)
        lh_y <- paste0("y",k+1-h)
        
        tmp_h <- tmp %>%
          mutate(hx = tmp[,h_x],
                 hy = tmp[,h_y],
                 lhx = tmp[,lh_x],
                 lhy = tmp[,lh_y])%>%
          mutate(xh = ifelse(x > xx, hx, lhx)) %>%
          mutate(yh = ifelse(x > xx, hy, lhy)) %>%
          mutate(xlh = ifelse(x > xx, lhx, hx)) %>%
          mutate(ylh = ifelse(x > xx, lhy, hy))%>%
          select(xh,yh,xlh,ylh)
        deto <- deto %>% select(-c(h_x,h_y,lh_x,lh_y))
        colnames(tmp_h) <- c(h_x,h_y,lh_x,lh_y)
        print(c(h_x,h_y,lh_x,lh_y))
        deto <- cbind(deto,tmp_h)
      }
    }
    nodus <- rbind(nodus,deto)
    k <- k + 1
    more <- more[!is.na(more[,tag_x]),] # reduce what left for the next loop
  }
  
  nodus%>%
    ggplot()+
    geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2))+
    geom_segment(aes(x = x2, xend = x3, y = y2, yend = y3))+
    geom_segment(aes(x = x3, xend = x4, y = y3, yend = y4))+
    geom_segment(aes(x = x4, xend = x5, y = y4, yend = y5))+
    geom_segment(aes(x = x5, xend = x6, y = y5, yend = y6))+
    # geom_segment(aes(x = x6, xend = x7, y = y6, yend = y7))+
    coord_fixed()
  
  # deux = trois = quatr = cinq <- NULL
  # deux <- nodes[is.na(nodes$x3),] # if the third point is missing, then they are a one segment wall
  # if(nrow(deux)> 0){
  #   deux <- deux %>%
  #     group_by(id_cell) %>%
  #     dplyr::mutate(xx = x2) %>%
  #     dplyr::mutate(yy = y2)
  #   deux <- deux %>%
  #     ungroup() %>%
  #     mutate(vertical = ifelse(x == xx, "true", "false")) %>%
  #     mutate(x1 = ifelse(x > xx, x, xx)) %>%
  #     mutate(x2 = ifelse(x > xx, xx, x)) %>%
  #     mutate(y1 = ifelse(x > xx, y, yy)) %>%
  #     mutate(y2 = ifelse(x > xx, yy, y)) %>%
  #     # Bug fix when wall is perfectly vertical
  #     mutate(y1 = ifelse(x == xx,
  #                        ifelse(y > yy, yy, y), y1)) %>%
  #     mutate(y2 = ifelse(x == xx,
  #                        ifelse(y > yy, y, yy), y2))
  # }
  # 
  # more <- nodes[!is.na(nodes$x3),] # a wall with at least more than one segment
  # if(nrow(more) > 0){
  #   trois <- more[is.na(more$x4),]
  #   trois <- trois %>%
  #     group_by(id_cell) %>%
  #     dplyr::mutate(xx = x3) %>%
  #     dplyr::mutate(yy = y3)
  #   trois <- trois %>%
  #     ungroup() %>%
  #     mutate(vertical = ifelse(x == xx, "true", "false")) %>%
  #     mutate(x1 = ifelse(x > xx, x, xx)) %>%
  #     mutate(x3 = ifelse(x > xx, xx, x)) %>%
  #     mutate(y1 = ifelse(x > xx, y, yy)) %>%
  #     mutate(y3 = ifelse(x > xx, yy, y)) %>%
  #     # Bug fix when wall is perfectly vertical
  #     mutate(y1 = ifelse(x == xx,
  #                        ifelse(y > yy, yy, y), y1)) %>%
  #     mutate(y3 = ifelse(x == xx,
  #                        ifelse(y > yy, y, yy), y3))
  #   more <- more[!is.na(more$x4),] # walls with at least more than two segments
  #   if(nrow(more) > 0){
  #     quatr <- more[is.na(more$x5),]
  #     
  #     quatr <- quatr %>%
  #       group_by(id_cell) %>%
  #       dplyr::mutate(xx = x4) %>%
  #       dplyr::mutate(yy = y4)
  #     quatr <- quatr %>%
  #       ungroup() %>%
  #       mutate(vertical = ifelse(x == xx, "true", "false")) %>%
  #       mutate(x1 = ifelse(x > xx, x, xx)) %>%
  #       mutate(x2 = ifelse(x > xx, x2, x3)) %>%
  #       mutate(y2 = ifelse(x > xx, y2, y3)) %>%
  #       mutate(x3 = ifelse(x > xx, x3, x2)) %>%
  #       mutate(y3 = ifelse(x > xx, y3, y2)) %>%
  #       mutate(x4 = ifelse(x > xx, xx, x)) %>%
  #       mutate(y1 = ifelse(x > xx, y, yy)) %>%
  #       mutate(y4 = ifelse(x > xx, yy, y)) %>%
  #       # Bug fix when wall is perfectly vertical
  #       mutate(y1 = ifelse(x == xx,
  #                          ifelse(y > yy, yy, y), y1)) %>%
  #       mutate(y4 = ifelse(x == xx,
  #                          ifelse(y > yy, y, yy), y4))
  #     more <- more[!is.na(more$x5),] # walls with at least more than three segments
  #     if(nrow(more) > 0){
  #       six <- more[is.na(more$x7),]
  #       
  #       six <- six %>%
  #         group_by(id_cell) %>%
  #         dplyr::mutate(xx = x6) %>%
  #         dplyr::mutate(yy = y6)
  #       six <- six %>%
  #         ungroup() %>%
  #         mutate(vertical = ifelse(x == xx, "true", "false")) %>%
  #         mutate(x1 = ifelse(x > xx, x, xx)) %>%
  #         mutate(x2 = ifelse(x > xx, x2, x4)) %>%
  #         mutate(y2 = ifelse(x > xx, y2, y4)) %>%
  #         mutate(x4 = ifelse(x > xx, x4, x2)) %>%
  #         mutate(y4 = ifelse(x > xx, y4, y2)) %>%
  #         mutate(x5 = ifelse(x > xx, xx, x)) %>%
  #         mutate(y1 = ifelse(x > xx, y, yy)) %>%
  #         mutate(y5 = ifelse(x > xx, yy, y)) %>%
  #         # Bug fix when wall is perfectly vertical
  #         mutate(y1 = ifelse(x == xx,
  #                            ifelse(y > yy, yy, y), y1)) %>%
  #         mutate(y4 = ifelse(x == xx,
  #                            ifelse(y > yy, y, yy), y5))
  #       more <- more[!is.na(more$x6),] # walls with at least more than three segments
  #       
  #       if(nrow(more)> 0){
  #         cinq <- more[is.na(more$x6),]
  #         
  #         cinq <- cinq %>%
  #           group_by(id_cell) %>%
  #           dplyr::mutate(xx = x5) %>%
  #           dplyr::mutate(yy = y5)
  #         cinq <- cinq %>%
  #           ungroup() %>%
  #           mutate(vertical = ifelse(x == xx, "true", "false")) %>%
  #           mutate(x1 = ifelse(x > xx, x, xx)) %>%
  #           mutate(x2 = ifelse(x > xx, x2, x4)) %>%
  #           mutate(y2 = ifelse(x > xx, y2, y4)) %>%
  #           mutate(x4 = ifelse(x > xx, x4, x2)) %>%
  #           mutate(y4 = ifelse(x > xx, y4, y2)) %>%
  #           mutate(x5 = ifelse(x > xx, xx, x)) %>%
  #           mutate(y1 = ifelse(x > xx, y, yy)) %>%
  #           mutate(y5 = ifelse(x > xx, yy, y)) %>%
  #           # Bug fix when wall is perfectly vertical
  #           mutate(y1 = ifelse(x == xx,
  #                              ifelse(y > yy, yy, y), y1)) %>%
  #           mutate(y4 = ifelse(x == xx,
  #                              ifelse(y > yy, y, yy), y5))
  #         more <- more[!is.na(more$x6),]
  #         
  #         
  #         
  #         stop(print("I was tired of coding the next segment ... please continue yourself"))}
  #     }
  #   }
  #   
  # }
  # 
  # nodes <- rbind(deux,trois,quatr,cinq)
  
  
  return(nodus)
}

write_vtk <- function(sim, path = "cross_section.vtk"){
  
  nodes <- sim$nodes
  nodes$id_point <- paste0(nodes$x,";",nodes$y)
  nodes$id_point1 <- paste0(nodes$x1,";",nodes$y1)
  nodes$id_point2 <- paste0(nodes$x2,";",nodes$y2)
  
  p <- nodes%>%
    filter(!duplicated(id_point))%>%
    mutate(z = 0)
  
  vtk <- '# vtk DataFile Version 2.0\n'
  vtk <- paste0(vtk, path,'\n')
  vtk <- paste0(vtk,'ASCII\n')
  vtk <- paste0(vtk,'DATASET POLYDATA\n')
  vtk <- paste0(vtk,'POINTS ',length(unique(nodes$id_point)),' float\n')
  point_data <- paste0(p$x,' ',p$y,' ',p$z, '\n', collapse = "")
  vtk <- paste0(vtk, point_data)
  
  vtk <- paste0(vtk, 'LINES ',nrow(nodes),' ',3*nrow(nodes),'\n')
  for(h in 1:nrow(nodes)){
    pos <- paste0('2 ',which(p$id_point == nodes$id_point1[h])-1,' ', which(p$id_point == nodes$id_point2[h])-1,'\n')
    vtk <- paste0(vtk,pos)
  }
  
  ppol <- NULL
  for(i in unique(nodes$id_cell)){
    tmp_cell <- nodes[nodes$id_cell == i,]
    
    n_coor <- nrow(tmp_cell)
    tmp_pol <- paste0(n_coor)
    for (j in 1:nrow(tmp_cell)) {
      tmp_pol <- paste0(tmp_pol,' ',which(p$id_point == tmp_cell$id_point[j])-1)
    }
    # print(tmp_pol)
    tmp_pol <- paste0(tmp_pol,'\n', collapse = "")
    ppol <- paste0(ppol, tmp_pol)
    }
  
    te <- strsplit(ppol,split = "[\n ]")
    de <- length(unlist(te))
  
  
  vtk <- paste0(vtk, 'POLYGONS ',length(unique(nodes$id_cell)),' ',de,'\n')
  vtk <- paste0(vtk,ppol)
  # vtk <- paste0(vtk, 'CELL_DATA ',length(unique(nodes$id_cell)),'\n',
  #               'SCALARS cell_scalars int 1\n',
  #               'LOOKUP_TABLE default\n')
  # ltd <- paste0(c(0:length(unique(nodes$id_cell))-1), '\n', collapse = "")
  # vtk <- paste0(vtk,ltd)

  cat(vtk, file= path)
  print(paste0(path, " has been created"))
}

aer_in_geom_xml <- function(sim, path = "~/Thesis/2020-02 GRANAR_3D/MECHA_GRANAR/Projects/GRANAR/in/Maize_Geometry.xml"){
  
  require(xml2)
  if (is.null(path)) {
    warning("No path specified")
  }
  if(is.null(sim$id_aerenchyma)){
    warning("No aerenchyma id specified")
  }else{
    id_aerenchyma <- sim$id_aerenchyma
  }
  xml <- read_xml(path)
  
  aer <- xml_children(xml_find_all(xml, "//aerenchyma_range"))
  
  # newbee <- 'aerenchyma id="0"'
  new_siblings <- paste0('aerenchyma id="',id_aerenchyma,'"')
  
  xml_add_sibling(aer, new_siblings)
  
  xml_remove(aer[1])
  
  path <- paste0(c(unlist(str_split(path, ".xml"))[1]),"_aer.xml")
  
  write_xml(xml, path)
  return(TRUE)
}


read_anatomy_xml <- function(path = "current_root.xml"){
  
  file <- read_xml(path)
  wallo <- xml_find_all(file, "//points")
  nodes <- NULL
  for (ch in xml_find_all(file, "//cell")) {
    id_cell <- as.numeric(unlist(xml_attrs(ch))[1])
    type <- as.numeric(unlist(xml_attrs(ch))[2])
    tmp_wall <- xml_children(xml_children(ch))
    tmp_wall <- as.numeric(unlist(xml_attrs(tmp_wall)))
    for (i in tmp_wall) {
      coor <- xml_children(wallo[[i+1]])
      for (k in 1:(length(coor)-1)) {
        x1 <- as.numeric(unlist(xml_attrs(coor[k])))[1]
        y1 <- as.numeric(unlist(xml_attrs(coor[k])))[2]
        
        x2 <- as.numeric(unlist(xml_attrs(coor[k+1])))[1]
        y2 <- as.numeric(unlist(xml_attrs(coor[k+1])))[2]
        
        tmp <- data.frame(x= x1,y = y1, x1 = x1, y1 = y1, x2 = x2, y2= y2, id_cell = id_cell, type = type, id_wall = i+1)
        nodes <- rbind(nodes,tmp)
      }
    }
  }
  nodes <- nodes%>%
    dplyr::group_by(id_cell) %>%
    dplyr::mutate(my = mean(y),
                  mx = mean(x),
                  atan = atan2(y-my, x - mx)) %>%
    dplyr::arrange(id_cell, atan)
  
  
  print(nodes%>%
    ggplot()+
    geom_segment(aes(x = x1, xend = x2, y = y1, yend= y2, colour = factor(type)))+
    coord_fixed())
  
  return(TRUE)
}


read_mecha_output <- function(dir = "C:/Users/heymansad/Documents/Thesis/2020-07 ICP_Aer/MECHA_GRANAR/Projects/GRANAR/out/M1v4/Root/Project_Test/results/"){
  
  K1 <- read.delim(paste0(dir,"Macro_prop_1,0.txt"))
  K2 <- read.delim(paste0(dir,"Macro_prop_2,1.txt"))
  K3 <- read.delim(paste0(dir,"Macro_prop_4,2.txt"))
  axial <- str_split(paste0(K1[4,1]), " ")[[1]]
  kx <- as.numeric(unlist(regmatches(axial,gregexpr("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?",axial))))[1]
  radial <- str_split(paste0(K1[5,1]), " ")[[1]]
  kr1 <- as.numeric(unlist(regmatches(radial,gregexpr("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?",radial))))[1]
  radial <- str_split(paste0(K2[5,1]), " ")[[1]]
  kr2 <- as.numeric(unlist(regmatches(radial,gregexpr("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?",radial))))[1]
  radial <- str_split(paste0(K3[5,1]), " ")[[1]]
  kr3 <- as.numeric(unlist(regmatches(radial,gregexpr("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?",radial))))[1]
  
  Ktot <- tibble(kx = kx,kr1 = kr1, kr2 = kr2, kr3 = kr3)
  return(Ktot)
}
  
  
  
  
  

