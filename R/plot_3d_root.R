
# Plot root in 3d
# rgl type

plot_granar3d <- function(root, col = "type"){
  if(col == "type"){
    name <- unique(root$type)
    for(n in 1:length(name)){
      root$type[root$type == name[n]] <- n
    }
    rgl.triangles(root$x1,root$x2,root$x3,col=root$type,alpha=0.5)
  }else{
    root$tmp <- root[,col]
    rgl.triangles(root$x1,root$x2,root$x3,col=root$tmp,alpha=0.5)
  }
}
