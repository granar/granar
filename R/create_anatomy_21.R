
#' @title GRANAR: Generate root anatomy in R
#'
#' This function creates a 2D root cross section anatomy based on anatomical features translated into the parameter of the model. It is the core of GRANAR.
#' @param path Path to the XML file containing the different parameters for the simulation. Not needed if 'parameter' is set.  Default = NULL.
#' @param parameters Table with the different parameters. Not needed if 'path' is set.  Default = NULL.
#' @param verbatim Display textual information aboutt he simulation. Default = NULL.
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

#source("~/Thesis/2020-02 GRANAR_3D/anatin_fun.R")

create_anatomy_2 <- function(path = NULL,  # PAth
                           parameters = NULL,
                           verbatim = F,
                           maturity_x = F,
                           paraview = F){

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

  t_1 <- Sys.time()
  t1 <- proc.time()

  # PARAMETERS -----
   if(verbatim) message("Loading parameters")
  # plant_type <- params$value[params$name == "planttype"]
   random_fact <- params$value[params$name == "randomness"] / 10 * params$value[params$name == "stele" & params$type == "cell_diameter"]
  # stele_diameter <- params$value[params$name == "stele" & params$type == "layer_diameter"]
  # n_aerenchyma_files <- params$value[params$name == "aerenchyma" & params$type == "n_files"]
   proportion_aerenchyma <- params$value[params$name == "aerenchyma" & params$type == "proportion"]
  # n_xylem_files <- params$value[params$name == "xylem" & params$type == "n_files"]
  # proto_meta_ratio <- params$value[params$name == "xylem" & params$type == "ratio"]
  # n_proto_xylem <- round(n_xylem_files*proto_meta_ratio)

  t2 <- proc.time()

  #//////////////////////////////////////////////////////////////////////////////////////////////////
  # INITIALIZE LAYERS -----

  if(verbatim) message("Creating cell layers")

  data_list <- cell_layer(params)
  layers <- data_list$layers
  all_layers <- data_list$all_layers

  center <- max(all_layers$radius)

  t3 <- proc.time()



  #//////////////////////////////////////////////////////////////////////////////////////////////////
  # CREATE CELLS ------

  if(verbatim) message("Creating cells")

  #all_layers$name[substr(all_layers$name, 1,6) == "cortex"] <- "cortex"
  all_cells <- create_cells(all_layers, random_fact)

  summary_cells <- ddply(all_cells, .(type), summarise, n_cells = length(angle))

  t4 <- proc.time()


  #//////////////////////////////////////////////////////////////////////////////////////////////////
  # CREATE XYLEM VESSELS -----
  # Create the xylem files
  # Get the extremes


  all_cells$type[substr(all_cells$type, 1,6) == "cortex"] <- "cortex"

  all_cells$id_group <- 0
  # ///////////////////////////////////////////////////////////////////////////////////


  if(verbatim) message("creating inter cellular spaces")
  # create inter-cellular space in cortex layer

  all_cells <- rondy_cortex(params, all_cells, center)


  all_cells %>%
    filter(id_group != 0)%>%
    ggplot()+
    geom_point(aes(x,y, colour = factor(id_group), shape = type), alpha = 0.5)+
    coord_fixed()+
    scale_shape_manual(values=1:nlevels(factor(all_cells$type)))+guides(colour = F)

  # Create vascular vessels for dicot or monocot
  if(verbatim) message("Creating xylem and phloem vessels")

  all_cells <- vascular(all_cells, params, layers, center)

  # all_cells$id_point <- paste0(round(all_cells$x,3),";",round(all_cells$y,3)) # 1 µm precision
  #
  # all_cells<- all_cells%>%
  #   # dplyr::group_by(id_point)%>%
  #   # dplyr::mutate(n = n())
  #   filter(!duplicated(id_point))%>%
  #   select(-id_point)

  # all_cells$type[all_cells$n != 1]

  all_cells %>%
    filter(id_group != 0 | type == "inter_cellular_space")%>%
    ggplot()+
    geom_point(aes(x,y, colour = factor(id_group), shape = type), alpha = 1)+
    coord_fixed()+
    scale_shape_manual(values=1:nlevels(factor(all_cells$type)))+guides(colour = F)

  t5 <- proc.time()

  #//////////////////////////////////////////////////////////////////////////////////////////////////
  # CREATE GEOMETRY ------
  if(verbatim) message("Creating the geometry")

  # Get the voronio data
  vtess <- deldir(all_cells$x, all_cells$y, digits = 8)
  if(is.null(vtess)){
    return(NULL)
    }

  # plot(vtess)
  # # Remove the ouside cells, to get the voronoi data straight
  # all_cells <- all_cells  %>%
  #   filter(type != "outside")

  if(verbatim) message("merging voronoi with cell type data")

  vorono_list <- cell_voro(all_cells, vtess, center)

  all_cells <- vorono_list$all_cells
  rs2 <- vorono_list$rs2

  t6 <- proc.time()

  rs1 <- rs2 %>%
    # dplyr::mutate(x = round(x,3),
    #               y = round(y,3),
    #               id_point = paste0(x,";",y))%>%
    dplyr::group_by(id_cell) %>%
    dplyr::mutate(my = mean(y),
                  mx = mean(x),
                  atan = atan2(y-my, x - mx)) %>%
    dplyr::arrange(id_cell, atan) #%>%
    #filter(!duplicated(atan))# %>% select(-id_point)

  # if(verbatim) message("endodermis structure")
  # endodermis rectangular shapes

  if(verbatim) message("Merging xylem vessels and cortex cells")
  rs1 <- smoothy_cells(rs1)

  rs1 <- rs1%>%
    dplyr::group_by(id_cell)%>%
    dplyr::filter(!duplicated(id_point))

  voiz <- rs1%>%
    mutate(id_point = paste0(round(x,3),";",round(y,3)))%>%
    dplyr::group_by(id_point)%>%
    dplyr::summarise(n = n())

  # message(paste0("a few possible mistake are possible around point", voiz$id_point[voiz$n < 2]))

  rs1 <- rs1[rs1$id_point %!in% voiz$id_point[voiz$n < 3] | rs1$type == "epidermis", ]

  if(verbatim) message("Merging inter cellular space")
  rs1 <- fuzze_inter(rs1)
  # rs1 <- rs1[rs1$type != "inter_cellular_space", ]

  if(verbatim){
  print(rs1%>%
          mutate(x = round(x,3),
                 y = round(y,3))%>%
          ggplot(aes(x, y, group=id_cell, fill = type))+
          geom_polygon(colour="white")+
          coord_fixed())
  }

  t7 <- proc.time()
  #//////////////////////////////////////////////////////////////////////////////////////////////////
  # CREATE AERENCHYMA -----
  if(verbatim){message("Killing cells to make aerenchyma")}

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

    if(verbatim){
      print(rs1%>%
              mutate(x = round(x,3),
                     y = round(y,3))%>%
              ggplot(aes(x, y, group=id_cell, fill = type))+
              geom_polygon(colour="white")+
              coord_fixed())
    }

    if(verbatim) message("simplify remaining cell walls between lacuna")
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

  # if(proportion_aerenchyma > 0){
  #   rs1$aer <- ifelse(rs1$id_cell %in% id_aerenchyma, "aer", "other")
  #   rs1$type[rs1$aer == "aer"] <- "aerenchyma"
  # }

  t8 <- proc.time()


  rs1 <- rs1[!is.na(rs1$x),]
  #//////////////////////////////////////////////////////////////////////////////////////////////////
  ## TIDY DATA ------

  if(verbatim) message("Tidying data before export")

  tt <- proc.time()
  # outputing the inputs
  output <- data.frame(io = "input", name = params$name, type = params$type, value = params$value)

  rs1$x <- round(rs1$x,3)
  rs1$y <- round(rs1$y,3)
  for (i in unique(rs1$id_cell)) {
    tmp <- rs1[rs1$id_cell == i,]
    tmp <- tmp[!is.na(tmp$x), ]
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

  one_cells <- rs1%>%
    filter(!duplicated(id_cell))# , !duplicated(type), !duplicated(id_group), !duplicated(area)
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
  TCA <- sum(all_cells$area.y[all_cells$order > 3])
  output <-rbind(output, data.frame(io="output", name="cortex_to_epidermis", type="layer_area",
                                    value = TCA))

  output <-rbind(output, data.frame(io="output", name="all", type="layer_area", value = sum(all_cells$area.y)))
  # output <-rbind(output, data.frame(io="output", name="aerenchyma", type="layer_area", value = (ini_cortex_area - cortex_area)))
  # output <-rbind(output, data.frame(io="output", name="aerenchyma", type="proportion", value = (ini_cortex_area - cortex_area)/ini_cortex_area))
  output <-rbind(output, data.frame(io="output", name="simulation", type="time", value = time))

  print(Sys.time()-t_1)

  rs1$sorting <- c(1:nrow(rs1))

  nodes <- vertex(rs1)

  # In the MECHA python script, to have unmature metaxylem vessels
  # Metaxylem elements are turned into stele cell type
  if(maturity_x){
    tmp_m <- mean(nodes$area[nodes$type == "cortex"])
    nodes$type[nodes$type == "xylem" & nodes$area > tmp_m ] <- "stele"
  }


  if(paraview){
    walls <- pv_ready(rs1)
    if(length(which(walls$keep_going == "...")) > 0){
      warning("walls should be described with more than 8 points")
    }
    wall_length <- walls%>%select((starts_with("x") | starts_with("y")) & ends_with(as.character(c(0:9))))%>%
      colnames()
    wally <- walls[!duplicated(walls[,wall_length]),] %>%
      dplyr::select(wall_length)
    wally$id_wall <- c(1:nrow(wally))
    walls <- merge(walls, wally, by= wall_length)
    walls <- walls %>%
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

  t9 <- proc.time()
  id_aerenchyma <- unique(nodes$id_cell[nodes$type %in% c("aerenchyma", "inter_cellular_space")])
  id_aerenchyma <- id_aerenchyma-1

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
              walls_nodes = walls,
              walls = wally,
              cells=all_cells,
              output = output,
              id_aerenchyma = id_aerenchyma))
}

write_anatomy_xml <- function(sim = NULL, path = NULL){

  if(is.null(sim)) warning("No simulation found. Please input a GRANAR simulation")
  if(is.null(path)) warning("No path found to save the XML file")

  if(length(sim$walls$x3) > 0){
    sim$nodes <- sim$walls_nodes
  }

  cellgroups <- data.frame(id_group = c(1, 2, 3, 3, 4, 5, 13, 16, 12, 11, 4, 4, 3),
                           type = c("exodermis", "epidermis", "endodermis", "passage_cell",  "cortex", "stele", "xylem", "pericycle", "companion_cell", "phloem", "inter_cellular_space", "aerenchyma", "cambium"))

  xml <- '<?xml version="1.0" encoding="utf-8"?>\n'
  xml <- paste0(xml, '<granardata>\n')

  # Write the Metadata
  xml <- paste0(xml, '\t<metadata>\n')
  xml <- paste0(xml, '\t\t<parameters>\n')
  xml <- paste0(xml,paste0('\t\t\t<parameter io="',sim$output$io,'" ',
                           'name="',sim$output$name,'" ',
                           'type="',sim$output$type,'" ',
                           'value="',sim$output$value,'"/>\n', collapse = ""))
  xml <- paste0(xml, '\t\t</parameters>\n')
  xml <- paste0(xml, '\t</metadata>\n')

  # Write the cells information
  xml <- paste0(xml, '\t<cells count="',nrow(sim$cells),'">\n')

  sim$nodes <- merge(sim$nodes, cellgroups, by="type")  %>%
    mutate(id_group = id_group.y)

  temp_wall <- ddply(sim$nodes, .(id_cell, id_group), summarise, walls = paste0('\t\t\t\t<wall id="',
                                                                                paste(id_wall-1, collapse='"/>\n\t\t\t\t<wall id="'),
                                                                                '"/>\n'))
  xml <- paste0(xml, paste0('\t\t<cell id="',temp_wall$id_cell-1, '" group="', temp_wall$id_group, '" truncated="false" >\n',
                            '\t\t\t<walls>\n', temp_wall$walls, '\t\t\t</walls>\n',
                            '\t\t</cell>\n', collapse=""))
  xml <- paste0(xml, '\t</cells>\n')


  # Write the walls information
  xml <- paste0(xml, '\t<walls count="',nrow(sim$walls),'">\n')
  walls <- sim$walls



  col_nam <- sim$walls%>%
    select((starts_with("x") | starts_with("y")) & ends_with(as.character(c(0:9))))%>%
    colnames()
  N <- max(parse_number(col_nam))
  begin <- tibble(tag1 = '\t\t<wall id="',
                  id_wall = sim$walls$id_wall-1,
                  tag2 = '" group="0" edgewall="false" >\n\t\t\t<points>\n')
  middle <- tibble(tag_x1 = '\t\t\t\t<point x="',
                   x1 = sim$walls$x1,
                   tag_y1 = '" y="',
                   y1 = sim$walls$y1,
                   tag_end1 = '"/>\n')
  for(k in 2:N){
    h <- k*2-1 # odd number
    tmp_coord <- sim$walls%>%
      select(all_of(col_nam[c(h,h+1)]))
    tmp_middle <- tibble(tag_x = '\t\t\t\t<point x="',
                         x = tmp_coord[,1],
                         tag_y = '" y="',
                         y = tmp_coord[,2],
                         tag_end = '"/>\n')
    tmp_col_name <- colnames(tmp_middle)
    colnames(tmp_middle) <- paste0(t(tmp_col_name), k)
    middle <- cbind(middle, tmp_middle)
  }
  taged_walls <- cbind(begin,middle)%>%
    mutate(tag_ending = '\t\t\t</points>\n\t\t</wall>\n')
  xml <- paste0(xml, paste0(t(taged_walls), collapse = ""))
  xml <- paste0(xml, '\t</walls>\n')
  xml <- str_remove_all(xml, '\t\t\t\t<point x=\"NA\" y=\"NA\"/>\n')

  print(cellgroups)
  xml <- paste0(xml, '\t<groups>\n')
  xml <- paste0(xml, '\t\t<cellgroups>\n')
  for(i in c(1:nrow(cellgroups))){
    xml <- paste0(xml, '\t\t\t<group id="',cellgroups$id_group[i],'" name="',cellgroups$type[i],'" />\n')
  }
  xml <- paste0(xml, '\t\t</cellgroups>\n')
  xml <- paste0(xml, '\t\t<wallgroups>\n')
  xml <- paste0(xml, '\t\t\t<group id="0" name="unassigned" />\n')
  xml <- paste0(xml, '\t\t</wallgroups>\n')
  xml <- paste0(xml, '\t</groups>\n')

  xml <- paste0(xml, '</granardata>')

  if(!is.null(path)){
    cat(xml, file = path)
    return(TRUE)
  }else{
    return(xml)
  }

}
