#' @title Generate a root cross-section
#'
#' Functions to generate root cross section anatomy, based on global parameters, such as the mean size of cells and number of cell layers.
#' @param path The path to the input file
#' @param parameters the input parameter table. Can be obtain by read_param_xml()
#' @param verbatim TRUE = Generate text to follow the simulation process/ FALSE = no text
#' @param maturity_x TRUE = meta-xylem are labeled as stele cell /FALSE = no change in cell labeling (default)
#' @param paraview TRUE = cell wall data is set to make 3D object/ FALSE = cell wall data is not compatible for 3D object
#' @keywords root
#' @import xml2
#' @import purrr
#' @import dplyr
#' @import tidyverse
#' @import deldir
#' @import sp
#' @import maptools
#' @import packcircles
#' @export
#' @examples
#' # Load input
#' params <- read_param_xml(path = system.file("extdata", "root_monocot.xml", package = "granar"))
#' # Generate anatomy
#' root = create_anatomy(parameters = params)
#' # Visualize the simulation output
#' plot_anatomy(root)
#' # Write the simulation output
#' write_anatomy_xml(sim = root, path = system.file("extdata", "current_root.xml", package = "granar"))
#'




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

  # set all cell center
  all_cells <- create_cells(all_layers, random_fact)
  # Get summary of cells
  summary_cells <- plyr::ddply(all_cells, plyr::.(type), summarise, n_cells = length(angle))
  # Label Cortex cells
  all_cells$type[grepl("cortex", all_cells$type)]<- "cortex"
  # Initialize id_group variable
  all_cells$id_group <- 0

  # Inclusion of the pith in the stele
  if(length(params$value[params$name == "pith" & params$type == "layer_diameter"]) > 0){
  all_cells <- make_pith(all_cells, params, center)
  }

  # Addition of intercellular space and reshape cortex layers
  # Sub optimal process, may take a while
  if(length(params$value[params$name =="inter_cellular_space"]) > 0){
    all_cells <- rondy_cortex(params, all_cells, center)
  }

  # Get the vascular system inside the stele
  # choose growth condition
  # if none --> do primary growth
  if(verbatim) message("Add vascular elements")
  if(length(params$value[params$name == "secondarygrowth"]) ==  0){
    all_cells <- vascular(all_cells, params, layers, center)
  } else if(params$value[params$name == "secondarygrowth"] == 0){
    all_cells <- vascular(all_cells, params, layers, center)
  } else if (params$value[params$name == "secondarygrowth"] == 1){
    # if sec growth, then do circle packing
    packing<-pack_xylem(all_cells, params, center)
    rm_stele <- all_cells%>%
      filter(type != "stele")
    new_cells <- rbind(packing, rm_stele)
    new_cells$id_cell <- 1:nrow(new_cells)
    all_cells <- new_cells
  }

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
    dplyr::arrange(id_cell, atan)%>%
    ungroup()

  rs1$id_point <- paste0(rs1$x,";",rs1$y)

  # Uniform cell by id_group
  if(verbatim) message("Smooth edge of large cells")
  rs1 <- smoothy_cells(rs1)

  if(verbatim) message("Merging inter cellular space")
  rs1 <- fuzze_inter(rs1)

  ini_cortex_area <- sum(all_cells$area[all_cells$type %in% c( "cortex" ,"exodermis" , # ,"endodermis"
                                                               "epidermis", "inter_cellular_space")])

  if(proportion_aerenchyma > 0){
    rs1 = clear_nodes(rs1)
    if(verbatim) message("remove cells for aerenchyma")
    rs1 <- aerenchyma(params, rs1)
  # simplify septa
    if(verbatim) message("simplify septa between aerenchyma lacuna")
    rs1 <- septa(rs1)
  }else{
    cortex_area <- ini_cortex_area
  }



  # hairy epidermis # add-on 27-02-2020
  #-----------------------------------------

  if(length(params$value[params$name == "hair"] )!= 0){
    if(params$value[params$name == "hair" & params$type == "n_files"] > 0 ){
      if(verbatim) message("Add root hair")
      rs1 <- root_hair(rs1, params, center)
    }
  }

  tt <- proc.time()
  # outputing the inputs
  output <- data.frame(io = "input", name = params$name, type = params$type, value = params$value)

  if(length(which(is.na(rs1$x)))>0 & verbatim){
    print("NA in cell coordinate ... ")
  }
  rs1 = clear_nodes(rs1)

  # Reset the ids of the cells to be continuous
  ids <- data.frame(id_cell = unique(rs1$id_cell))
  ids$new <- c(1:nrow(ids))
  rs1 <- merge(rs1, ids, by="id_cell")
  rs1$id_cell <- rs1$new

  if(proportion_aerenchyma > 0){
    if(verbatim) message("create id_aerenchyma vector")
    id_aerenchyma <- unique(rs1$id_cell[rs1$aer == "aer"])
  }else{id_aerenchyma <- NA}

  all_cells <- merge(all_cells, ids, by="id_cell")
  all_cells$id_cell <- all_cells$new

  mX <- mean(rs1$area[rs1$type == "xylem"])
  if(params$value[params$name == "planttype"] == 1){
    if(verbatim) message("for monocot, if xylem is above average, it is labeled as metaxylem")
    rs1$type[rs1$type == "xylem" & rs1$area > mX] <- "metaxylem"
  }
  one_cells <- rs1%>%
    filter(!duplicated(id_cell))# , !duplicated(type), !duplicated(id_group), !duplicated(area)

  all_cells <- merge(all_cells, one_cells, by = "id_cell")

  # adding the outputs by cell layers
  out <- plyr::ddply(all_cells, plyr::.(type.y), summarise, n_cells=length(type.y),
               layer_area = sum(area.y),
               cell_area = mean(area.y)) %>%
    mutate(name = type.y) %>%
    dplyr::select(-type.y) %>%
    tidyr::gather(key = "type", value = "value", n_cells, layer_area, cell_area) %>%
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

 # comment
  if(paraview){
    walls <- pv_ready(rs1)
    wall_length <- walls%>%select(-x, -y, -xx, -yy)%>% # ends_with(as.character(c(0:9)))
      select(starts_with("x"), starts_with("y"))%>%
      colnames()
    if(verbatim){
      print(wall_length)
    }

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

`%!in%` <- compose(`!`, `%in%`)
