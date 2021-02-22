#' Plot root anatomy
#'
#' This function plot the results of a 2D root cross section anatomy simulation
#' @param sim the simulation objkect, returned by 'create_anatomy.R'
#' @param col Parameter to choose for the coloring of the cells. Accepted arguments are 'type', 'area', 'dist', 'id_cell', "segment' and 'angle'. Default = 'type'
#' @param leg Display the legend; Default= TRUE
#' @param apo_bar Display apolastic barrier when col = "segment". 1 endodermal casparian strip, 2 fully suberized endodermis, 3 fully suberized endodermis and an exodermal casparian strip, and 4 exodermis and endodermis are fully suberized.
#' @keywords root
#' @export
#' @import viridis
#' @examples
#' create_anatomy("PATH_TO_XLM_FILE")
#'


plot_anatomy <- function(sim=NULL,
                         col = "type",
                         leg = T,
                         apo_bar = 2){

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

    if(apo_bar != 0){
    if(apo_bar %in% c(2,4)){
      if(apo_bar == 4){
        pl <- pl + geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), size = 1, colour = "red",
                                data = sim$nodes%>%
                                  filter(type == "exodermis"))

      }
      pl <- pl + geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), size = 1, colour = "red",
                              data = sim$nodes%>%
                                filter(type == "endodermis"))
    }
    if(apo_bar %in% c(1,3)){
      nod_apo <- sim$nodes
      cx <- unique(nod_apo$mx[nod_apo$id_cell == 1])
      cy <- unique(nod_apo$my[nod_apo$id_cell == 1])
      nod_apo <- nod_apo%>%
        mutate(slo = atan((y2-y1)/(x2-x1)),
               SLO = atan((y1-cy)/(x1-cy)),
               d = (slo-SLO))

      if(apo_bar == 1){
      pl <- pl + geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), size = 1.2, colour = "red", data = nod_apo%>%
                                filter(type == "endodermis"
                                       ,d < 0.4 & d > - 0.4))
      }else{
        pl <- pl + geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), size = 1, colour = "red",
                                data = sim$nodes%>%
                                  filter(type == "endodermis"))+
                  geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), size = 1.2, colour = "red", data = nod_apo%>%
                                filter(type == "exodermis"
                                       ,d < 0.4 & d > - 0.4))
      }
    }
    }
  }

  return(pl)
}

