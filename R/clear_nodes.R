

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
      pol <- Polygon(tmp[, c("x","y")])
      rs1$area[rs1$id_cell == i] <-  pol@area
      if(pol@area == 0){
         rs1 <- rs1[rs1$id_cell != i,]
      }
    }
  }

  return(rs1)
}
