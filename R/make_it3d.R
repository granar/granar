# GRANAR 3D
# feb 2020
# Heymans Adrien

make_it3D <- function(sim, n_cell){
require(geometry)

  col_vect <- c("lightgoldenrod", "darkgoldenrod")
  col_vect <- rep(col_vect, n_cell)
  root <- NULL
  k <- 0
  for (i in unique(sim$nodes$id_cell)) {
    is.xylem <- sim$nodes$type[sim$nodes$id_cell == i][1] == "xylem"
    k <- k + 1
    te <- data.frame(x = sim$nodes$x[sim$nodes$id_cell == i], y = sim$nodes$y[sim$nodes$id_cell == i] )
    te$z  <- runif(1,-0.02,0.05)
    len <- te$z+0.1+runif(1,-0.05,0.05)
    if(is.xylem){
      te$z  <- runif(1,-0.1,-0.05)
      len <- te$z+1+runif(1,-0.1,0.1)
    }
    te <- rbind(te, te%>%mutate(z = len))
    te <- as.matrix(te)
    A.surf <- t(convhulln(te))
    A.vol <- convhulln(te, output.options = "FA")$vol
    tmp <- tibble(x1 = te[A.surf,1], x2 = te[A.surf,2], x3 = te[A.surf,3], col=col_vect[1],
                  id_cross = i, axis = 1, ID = k, volume = A.vol, type = sim$nodes$type[sim$nodes$id_cell==i][1],
                  area_cross = sim$nodes$area[sim$nodes$id_cell == i][1])
    root <- rbind(root, tmp)
    if(is.xylem){next}
    for(l in 2:n_cell){
      k <- k + 1
      te <- data.frame(x = sim$nodes$x[sim$nodes$id_cell == i], y = sim$nodes$y[sim$nodes$id_cell == i] )
      te$z  <- unique(max(tmp$x3))
      len <- te$z+0.1+runif(1,-0.05,0.05)
      if(sim$nodes$type[sim$nodes$id_cell == i][1] == "xylem"){
        len <- te$z+1+runif(1,-0.05,0.05)
      }
      te <- rbind(te, te%>%mutate(z = len))
      te <- as.matrix(te)
      A.surf <- t(convhulln(te))
      A.vol <- convhulln(te, output.options = "FA")$vol
      tmp <- tibble(x1 = te[A.surf,1], x2 = te[A.surf,2], x3 = te[A.surf,3], col=col_vect[l],
                    id_cross = i, axis = l, ID = k, volume = A.vol, type = sim$nodes$type[sim$nodes$id_cell==i][1],
                    area_cross = sim$nodes$area[sim$nodes$id_cell == i][1])
      root <- rbind(root, tmp)
    }
  }
  return(root)

}
