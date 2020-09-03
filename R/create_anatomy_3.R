
library(tidyverse)
library(plyr)
library(xml2)
library(deldir)
library(sp)
source("./read_param_xml.R")
source("./io_function.R")

params <- read_param_xml("~/Thesis/2020-07 ICP_Aer/Zea_mays_2020.xml")
# Let's try again

create_anatomy_3 <- function(path = NULL,  # path to xml file
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

  # set initial time
  t1 <- proc.time()
  # set the random factor
  random_fact <- params$value[params$name == "randomness"] / 10 * params$value[params$name == "stele" & params$type == "cell_diameter"]

  data_list <- cell_layer(params)
  layers <- data_list$layers # layers: cell_type, diameter, n_layer, order
  all_layers <- data_list$all_layers # expand layers
  center <- max(all_layers$radius) # center of the cross section

  # layer time
  t2 <- proc.time()

  #set all cell center
  all_cells <- create_cells(all_layers, random_fact)
  summary_cells <- ddply(all_cells, .(type), summarise, n_cells = length(angle))
  all_cells$type[substr(all_cells$type, 1,6) == "cortex"] <- "cortex"
  all_cells$id_group <- 0

  # get the vascular system inside the stele
  all_cells <- vascular(all_cells, params, layers, center)

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

  # rs1%>%
  #   ggplot()+
  #   geom_polygon(aes(x,y, group = id_cell, fill = type ), colour = "white")+
  #   coord_fixed()+
  #   theme_classic()

  rs1$id_point <- paste0(rs1$x,";",rs1$y)

  # ! The function that mess the code !
  rs1 <- smoothy_cells(rs1)

  neib <- rs1$id_point[rs1$id_cell == 59]
  nebo <- rs1$id_cell[rs1$id_cell != 59 & rs1$id_point %in% neib]

  tmp_rs <- rs1
  tmp_rs$type[tmp_rs$id_cell == 59] <- "one"
  tmp_rs$type[tmp_rs$id_cell %in% nebo] <- "neibo"

  tmp_rs%>%
    ggplot()+
    geom_polygon(aes(x,y, group = id_cell, fill = type ), colour = "white")+
    coord_fixed()+
    theme_classic()

}
