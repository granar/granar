

#' @title Create node dataframe of null cells
#'
#' @param rs1 node dataframe
#' @keywords root
#' @export
#'

clear_nodes <- function(rs1){

  rs1 <- rs1[!is.na(rs1$x), ]
  for (i in unique(rs1$id_cell)) {
    tmp <- rs1[rs1$id_cell == i,]
    if(nrow(tmp)> 0){
      subset_points <- rbind(tmp[, c("x","y")],tmp[1, c("x","y")])
      poly_matrix <- matrix(c(subset_points$x, subset_points$y), ncol = 2)
      # Create an sf object for the subset points
      subset_sf <- sf::st_polygon(list(poly_matrix))
      rs1$area[rs1$id_cell == i] <-  sf::st_area(subset_sf)
      if(sf::st_area(subset_sf) == 0){
        rs1 <- rs1[rs1$id_cell != i,]
      }
    }
  }

  return(rs1)
}
