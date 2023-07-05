#' @title Smooth the edge of cells
#'
#' @param data node data frame
#' @keywords root
#' @export
#'

smoothy_cells <- function(data){

  for(i in c(1:max(data$id_group))){
    temp <- data[data$id_group == i,]
    data <- data[data$id_group != i,]
    if(nrow(temp)> 0){
      temp <- concavety(temp)
      data <- rbind(data,temp)%>%dplyr::as_tibble()
    }
  }

  data <- data%>%
    dplyr::group_by(id_cell)%>%
    dplyr::filter(!duplicated(id_point))%>%
    ungroup()
  voiz <- data%>%
    dplyr::group_by(id_point)%>%
    dplyr::summarise(n = n())%>%
    ungroup()

  data <- data[data$id_point %!in% voiz$id_point[voiz$n < 3] | data$type == "epidermis", ]

  return(data)
}
