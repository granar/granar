#' @title Make a 3D plot (rgl) root of a few cells
#'
#'
#' @param root The 3d root data
#' @param col what should be highlited
#' @keywords root 3D
#' @export
#'

plot_granar3d <- function(root, col = "type"){
  if(col == "type"){
    name <- unique(root$type)
    for(n in 1:length(name)){
      root$type[root$type == name[n]] <- n
    }
    rgl::rgl.triangles(root$x1,root$x2,root$x3,col=root$type,alpha=0.5)
  }else{
    root$tmp <- root[,col]
    rgl::rgl.triangles(root$x1,root$x2,root$x3,col=root$tmp,alpha=0.5)
  }
}
