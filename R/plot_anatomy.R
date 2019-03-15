#' Plot root anatomy
#'
#' This function plot the results of a 2D root cross section anatomy simulation
#' @param sim the simulation objkect, returned by 'create_anatomy.R'
#' @param col Parameter to choose for the coloring of the cells. Accepted arguments are 'type', 'area', 'dist', 'id_cell' and 'angle'. Default = 'type'
#' @param leg Display the legend; Default= TRUE
#' @keywords root
#' @export
#' @import viridis
#' @examples
#' create_anatomy("PATH_TO_XLM_FILE")
#'


plot_anatomy <- function(sim=NULL,
                         col = "type",
                         leg = T){

  pl <- ggplot(sim$nodes) +
    geom_polygon(aes_string("x", "y", group="id_cell", fill=col), colour="white") +
    theme_classic() +
    coord_fixed() +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())

  if(!col %in% c("type", "cell_group")){
    pl <- pl + scale_fill_viridis()
  }

  if(!leg){
    pl <- pl + theme(legend.position="none")
  }


  return(pl)

}
