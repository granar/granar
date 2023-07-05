#' @title Make Septa in between aerenchyma lacuna
#'
#'
#' @param rs1 node dataframe
#' @keywords root
#' @export
#' @examples
#' # rs <- septa(rs1)
#'
#'

septa <- function(rs1){
  data <- rs1
  data$id_point <- paste0(data$x,";",data$y)
  data$neib1 = data$neib2 <- NA
  data$neib1_type = data$neib2_type <- NA

  for(i in which(data$type == "aerenchyma")){
    ne <- unique(data$id_cell[data$id_point == data$id_point[i] & data$id_cell != data$id_cell[i]])
    data$neib1[i] <- ne[1]
    data$neib1_type[i] <- data$type[data$id_cell == ne[1]][1]
    data$neib2[i] <- ne[2]
    data$neib2_type[i] <- data$type[data$id_cell == ne[2]][1]
    if(!is.na(ne[3])){message("quadri point")
      print(ne)}
  }
  # all triple points and all points next to cortex cells
  must <- data[!is.na(data$neib2),]
  must <- rbind(must, data[data$neib1_type != "aerenchyma",])

  # add some noise
  noise <- data[!(data$neib1_type != "aerenchyma" | !is.na(data$neib2)),]
  noise_point <- unique(noise$id_point)
  n_noise <- round(length(noise_point)*0.05)
  noise_keep <- sample(noise_point, n_noise, replace=F)
  noise <- noise[noise$id_point %in% noise_keep,]
  must <- rbind(must, noise)

  data <- data[data$type != "aerenchyma",]
  data <- rbind(data,must)%>%
    arrange(id_cell,atan)

  return(data%>%select(colnames(rs1)))
}
