
#' write a .obj from the generated root cross section
#'
#' @param sim List generated from create_anatomy()
#' @param path Where to write the generated export
#' @param membrane create cell membrane if membrane = T, create cell wall if membrane = F
#' @keywords root anatomy
#' @export
#' @examples
#' write_CT_OBJ(sim, path = "cross_section.obj)
#'

write_CT_OBJ <- function(sim, path = "cross_section.obj", membrane = T){

  nodes <- sim$nodes

  if(membrane){
  nodes$cx <- (nodes$x + nodes$mx)/2
  nodes$cy <- (nodes$y + nodes$my)/2

  for(i in 1:4){
    nodes$cx <- (nodes$x + nodes$cx)/2
    nodes$cy <- (nodes$y + nodes$cy)/2
  }
  nodes$x <- nodes$cx
  nodes$y <- nodes$cy
  }

  nodes$id_point <- paste0(nodes$x,";",nodes$y)

  p <- nodes%>%
    filter(!duplicated(id_point))%>%
    mutate(z = 0)

  obj <- '# Simple Wavefront file\n'
  point_data <- paste0('v ',p$x,' ',p$y,' ',p$z, '\n', collapse = "")
  obj <- paste0(obj, point_data)

  if(membrane){
    ppol <- NULL
    for(i in unique(nodes$id_cell)){
      tmp_cell <- nodes[nodes$id_cell == i,]
      # create a cell face
      tmp_pol <- 'f'
      for (j in 1:nrow(tmp_cell)) {
        tmp_pol <- paste0(tmp_pol," ", which(p$id_point == tmp_cell$id_point[j]))
      }
      # print(tmp_pol)
      tmp_pol <- paste0(tmp_pol,'\n', collapse = " ")
      ppol <- paste0(ppol, tmp_pol)
    }
  }else{
    nodes$id2_points <- paste0(nodes$x2,";",nodes$y2)
    nodes$id1_points <- paste0(nodes$x1,";",nodes$y1)

    ppol <- paste0('l ',match(nodes$id1_points, p$id_point),' ',
                   match(nodes$id2_points, p$id_point), '\n', collapse = "")
  }
  obj <- paste0(obj,ppol)

  cat(obj, file= path)
  print(paste0(path, " has been created"))
}
