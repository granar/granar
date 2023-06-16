#' @title Place one line of cell center
#'
#' Create a list
#' @param params The input dataframe
#' @keywords root layer
#' @export
#' @examples
#' data_list <- cell_layer(params)
#'

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
