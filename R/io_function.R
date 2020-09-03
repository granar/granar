

# io_function for create_anatomy_3


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

