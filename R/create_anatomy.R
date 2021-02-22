

create_anatomy <- function(path = NULL,  # path to xml file
                           parameters = NULL,
                           verbatim = F,
                           maturity_x = F,
                           paraview = T){
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

  # set initial time
  t_1 <- Sys.time()
  t1 <- proc.time()
  # set the random factor
  random_fact <- params$value[params$name == "randomness"] / 10 * params$value[params$name == "stele" & params$type == "cell_diameter"]
  proportion_aerenchyma <- params$value[params$name == "aerenchyma" & params$type == "proportion"]

  data_list <- cell_layer(params)
  layers <- data_list$layers # layers: cell_type, diameter, n_layer, order
  all_layers <- data_list$all_layers # expand layers
  center <- max(all_layers$radius) # center of the cross section

  # layer time
  t2 <- proc.time()

  #set all cell center
  all_cells <- create_cells(all_layers, random_fact)
  summary_cells <- plyr::ddply(all_cells, .(type), summarise, n_cells = length(angle))
  all_cells$type[substr(all_cells$type, 1,6) == "cortex"] <- "cortex"

  all_cells$id_group <- 0

  if(length(params$value[params$name == "pith" & params$type == "layer_diameter"]) > 0){
  all_cells <- make_pith(all_cells, params, center)
  }

  # get the vascular system inside the stele
  if(length(params$value[params$name =="inter_cellular_space"]) > 0){
    all_cells <- rondy_cortex(params, all_cells, center)
  }
  all_cells <- vascular(all_cells, params, layers, center)


  all_cells%>%
    ggplot()+
    geom_point(aes(x,y, colour = type), alpha = 0.3)+
    coord_fixed()+
    theme_classic()

  # Get the voronio data
  vtess <- deldir(all_cells$x, all_cells$y, digits = 8)
  if(is.null(vtess)){return(NULL)}
  vorono_list <- cell_voro(all_cells, vtess, center)
  all_cells <- vorono_list$all_cells
  rs2 <- vorono_list$rs2

  rs1 <- rs2 %>%
    dplyr::group_by(id_cell) %>%
    dplyr::mutate(my = mean(y),
                  mx = mean(x),
                  atan = atan2(y-my, x - mx)) %>%
    dplyr::arrange(id_cell, atan)

  rs1$id_point <- paste0(rs1$x,";",rs1$y)

  rs1 <- rs1%>%
    dplyr::group_by(id_cell)%>%
    dplyr::filter(!duplicated(id_point))%>%
    ungroup()

  rs1%>%
    ggplot()+
    geom_polygon(aes(x,y, group = id_cell, fill = type), colour = "white")+
    coord_fixed()+
    theme_classic()

  rs1 <- smoothy_cells(rs1)

  # message(paste0("a few possible mistake are possible around point", voiz$id_point[voiz$n < 2]))

  if(verbatim) message("Merging inter cellular space")
  rs1 <- fuzze_inter(rs1)
  # rs1 <- rs1[rs1$type != "inter_cellular_space", ]

  ini_cortex_area <- sum(all_cells$area[all_cells$type %in% c( "cortex" ,"exodermis" , # ,"endodermis"
                                                               "epidermis", "inter_cellular_space")])
  saved_rs1 <- rs1
  rs1 <- saved_rs1
  if(proportion_aerenchyma > 0){

    for (i in unique(rs1$id_cell)) {
      tmp <- rs1[rs1$id_cell == i,]
      tmp <- tmp[!is.na(tmp$x), ]
      pol <- Polygon(tmp[, c("x","y")])
      rs1$area[rs1$id_cell == i] <-  pol@area
    }
    # make aerenchyma
    rs1 <- aerenchyma(params, rs1)

    # simplify septa
    rs1 <- septa(rs1)
    # id_aerenchyma <- unique(septum$id_cell)
  }else{cortex_area <- ini_cortex_area}


  # hairy epidermis # add-on 27-02-2020
  #-----------------------------------------

  if(length(params$value[params$name == "hair"] )!= 0){
    if(params$value[params$name == "hair" & params$type == "n_files"] > 0 ){
      rs1 <- root_hair(rs1, params, center)
    }
  }

  tt <- proc.time()
  # outputing the inputs
  output <- data.frame(io = "input", name = params$name, type = params$type, value = params$value)

  if(length(which(is.na(rs1$x)))>0){
    print("NA in cell coordinate ... ")
  }
  rs1 <- rs1[!is.na(rs1$x), ]
  for (i in unique(rs1$id_cell)) {
    tmp <- rs1[rs1$id_cell == i,]
    if(nrow(tmp)> 0){
      pol <- Polygon(tmp[, c("x","y")])
      rs1$area[rs1$id_cell == i] <-  pol@area
      if(pol@area == 0){
        print("cell_area = 0")
        print("this id_cell will be removed")
        rs1 <- rs1[rs1$id_cell != i,]
      }
    }
  }

  # Reset the ids of the cells to be continuous
  ids <- data.frame(id_cell = unique(rs1$id_cell))
  ids$new <- c(1:nrow(ids))
  rs1 <- merge(rs1, ids, by="id_cell")
  rs1$id_cell <- rs1$new

  if(proportion_aerenchyma > 0){
    id_aerenchyma <- unique(rs1$id_cell[rs1$aer == "aer"])
  }else{id_aerenchyma <- NA}

  all_cells <- merge(all_cells, ids, by="id_cell")
  all_cells$id_cell <- all_cells$new

  mX <- mean(rs1$area[rs1$type == "xylem"])
  if(params$value[params$name == "planttype"] == 1){
    rs1$type[rs1$type == "xylem" & rs1$area > mX] <- "metaxylem"
  }
  one_cells <- rs1%>%
    filter(!duplicated(id_cell))# , !duplicated(type), !duplicated(id_group), !duplicated(area)



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
  TCA <- sum(all_cells$area.y[all_cells$order > 3])
  output <-rbind(output, data.frame(io="output", name="cortex_alive_to_epidermis", type="layer_area",
                                    value = TCA))

  output <-rbind(output, data.frame(io="output", name="all", type="layer_area", value = sum(all_cells$area.y)))
  # output <-rbind(output, data.frame(io="output", name="aerenchyma", type="layer_area", value = (ini_cortex_area - cortex_area)))
  # output <-rbind(output, data.frame(io="output", name="aerenchyma", type="proportion", value = (ini_cortex_area - cortex_area)/ini_cortex_area))
  output <-rbind(output, data.frame(io="output", name="simulation", type="time", value = time))



  rs1$sorting <- c(1:nrow(rs1))

  nodes <- vertex(rs1)
  nodes <- nodes[!is.na(nodes$x), ]

  # In the MECHA python script, to have unmature metaxylem vessels
  # Metaxylem elements are turned into stele cell type
  if(maturity_x){
    tmp_m <- mean(nodes$area[nodes$type == "cortex"])
    nodes$type[nodes$type == "metaxylem" & nodes$area > tmp_m ] <- "stele"
  }


  if(paraview){
    walls <- pv_ready(rs1)
    wall_length <- walls%>%select(-x, -y, -xx, -yy)%>% # ends_with(as.character(c(0:9)))
      select(starts_with("x"), starts_with("y"))%>%
      colnames()
    print(wall_length)
    wally <- walls[!duplicated(walls[,wall_length]),] %>%
      dplyr::select(wall_length)
    wally$id_wall <- c(1:nrow(wally))
    walls <- merge(walls, wally, by= wall_length)
    walls <- walls %>%
      # filter(!duplicated(id_wall))%>% # 30/09/2020
      arrange(sorting)

  }else{
    wally <- nodes[!duplicated(nodes[,c('x1', 'x2', 'y1', 'y2')]),] %>%
      dplyr::select(c(x1, x2, y1, y2))

    wally$id_wall <- c(1:nrow(wally))
    walls <- wally

    nodes <- merge(nodes, walls, by=c("x1", "x2", "y1", "y2"))
    nodes <- nodes %>%
      arrange(sorting)
  }

  id_aerenchyma <- unique(nodes$id_cell[nodes$type %in% c("aerenchyma", "inter_cellular_space")])
  id_aerenchyma <- id_aerenchyma-1

  print(Sys.time()-t_1)

  return(list(nodes = nodes,
              walls_nodes = walls,
              walls = wally,
              cells=all_cells,
              output = output,
              id_aerenchyma = id_aerenchyma))

}
