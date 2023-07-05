
#' @title Compute the sum of the cell areas
#'
#'
#' @param cells the dataframe containing the cells to sum up
#' @keywords root
#' @export
#'

sum_area <- function(cells){
  area_tmp <- cells %>%
    dplyr::group_by(id_cell)%>%
    dplyr::summarise(area = mean(area))%>%
    ungroup()

  area_tmp <- area_tmp%>%
    plyr::arrange(area, decreasing = T)
  are <- sum(area_tmp$area)
  return(are)

}
