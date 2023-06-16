#' @title Get the cooridnates for every nodes in the voronoi
#'
#'
#' @param all_cells The cellular dataframe
#' @param vtess the voronoi data
#' @param center The cross-section center
#' @keywords root
#' @export
#' @examples
#'
#'

cell_voro <- function(all_cells, vtess, center){

  # Get the size of the cells
  ids <- all_cells$id_cell
  all_cells$area <- NA
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
