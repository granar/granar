#' Plot root anatomy
#'
#' This function plot the results of a 2D root cross section anatomy simulation
#' @param sim the simulation objkect, returned by 'create_anatomy.R'
#' @param col Parameter to choose for the coloring of the cells. Accepted arguments are 'type', 'area', 'dist', 'id_cell', "segment' and 'angle'. Default = 'type'
#' @param leg Display the legend; Default= TRUE
#' @param apo_bar Display apolastic barrier when col = "segment". 1 endodermal casparian strip, 2 fully suberized endodermis, 3 fully suberized endodermis and an exodermal casparian strip, and 4 exodermis and endodermis are fully suberized.
#' @param phi_thck Display phi thickening on the defined number of layer after the endodermis
#' @keywords root
#' @export
#' @import viridis
#' @import ggplot2
#' @examples
#' # sim = create_anatomy(path = "PATH_TO_XLM_FILE")
#' # plot_anatomy(sim, col = "type")
#' # plot_anatomy(sim, col = "segment", apo_bar = 3)
#'


plot_anatomy <- function(sim=NULL,
                         col = "type",
                         leg = T,
                         apo_bar = 0,
                         phi_thck = 0){

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
    pl <- pl + viridis::scale_fill_viridis()
  }

  if(!leg){
    pl <- pl + theme(legend.position="none")
  }
  if(col == "segment"){
    pl <- ggplot()+
      geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), data = sim$nodes)+
      theme_classic()+
      coord_fixed()+
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())
    # Addition of the hydrophobic barriers

    nod_apo <- sim$nodes
    cx <- unique(nod_apo$mx[nod_apo$id_cell == 1])
    cy <- unique(nod_apo$my[nod_apo$id_cell == 1])
    nod_apo <- nod_apo%>%
      mutate(slo = atan((y2-y1)/(x2-x1)),
             SLO = atan((y1-cy)/(x1-cy)),
             d = (slo-SLO),
             x_mid = (x1+x2)/2,
             y_mid = (y1+y2)/2)
    if(apo_bar == 0){}else{

      if(apo_bar == 1){
        pl <- pl + geom_point(aes(x = x_mid, y = y_mid), size = 1.2, colour = "red", data = nod_apo%>%
                                  filter(type == "endodermis"
                                         ,d < 0.4 & d > - 0.4| d > 2.8 | d < -2.8))
      }
      if(apo_bar >= 2){
        pl <- pl + geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), size = 1, colour = "red",
                                data = sim$nodes%>%
                                  filter(type == "endodermis"))
      }
      if(apo_bar == 3){
        pl = pl + geom_point(aes(x = x_mid, y = y_mid), size = 1.2, colour = "red", data = nod_apo%>%
                       filter(type == "exodermis"
                              ,d < 0.4 & d > - 0.4| d > 2.8 | d < -2.8))
      }
      if(apo_bar == 4){
        pl <- pl + geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), size = 1, colour = "red",
                                data = sim$nodes%>%
                                  filter(type == "exodermis"))

      }
    }

    if(phi_thck == 0){}else{
      id_endo_layer = unique(nod_apo$id_layer[nod_apo$type == "endodermis"])
      pl <- pl + geom_point(aes(x = x_mid, y = y_mid), size = 1.2, colour = "red", data = nod_apo%>%
                             filter(id_layer > id_endo_layer,
                                    id_layer <= id_endo_layer+phi_thck,
                                    type != "inter_cellular_space",
                                    type != "aerenchyma",
                                    d < 0.4 & d > - 0.4 | d > 2.8 | d < -2.8))
    }
  }



  return(pl)
}

