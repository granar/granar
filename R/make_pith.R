#' @title Make pith inside the stele
#'
#' remove stele cells and place pith cells
#' @param all_cells The cellular dataframe
#' @param params The input dataframe
#' @param center The cross section center location
#' @keywords root
#' @export
#' @examples
#' all_cells <- make_pith(all_cells, params, center)
#'


make_pith <- function(all_cells, params, center){
  if(params$value[params$name == "pith"][1] > 0){
    pith_size <- params$value[params$name == "pith" & params$type == "layer_diameter"]/2
    pcell <- params$value[params$name == "pith" & params$type == "cell_diameter"]
  }else{pith_size <- 0}

  if(pith_size > 0){


    xylem <- all_cells%>%
      filter(type == "xylem")

    xylem <- xylem%>%
      dplyr::group_by(id_group)%>%
      dplyr::mutate(mx = mean(x),
                    my = mean(y),
                    euc = sqrt((mx-center)^2+(my - center)^2))
    inner <- unique(xylem$id_group[xylem$euc < pith_size])
    all_cells <- all_cells[all_cells$type != "xylem" | all_cells$id_group %!in% inner,]


    n_pith_lay <- round(1+(pith_size-pcell/2)/pcell)
    pith_layer <- data.frame(name="stele",
                             n_layers=rep(n_pith_lay, n_pith_lay),
                             cell_diameter=pcell,
                             order = 0.5)

    pith_layer <- layer_info(pith_layer)
    new_cells <- create_cells(all_layers = pith_layer, random_fact = 0.001)
    new_cells%>%
      ggplot()+geom_point(aes(x,y))+
      coord_fixed()
    new_center <- mean(new_cells$x[new_cells$angle == 0], new_cells$y[new_cells$angle == 0])

    new_cells$x <- new_cells$x-new_center+center
    new_cells$y <- new_cells$y-new_center+center
    new_cells$id_group <- 0

    all_cells <- all_cells[sqrt((all_cells$x-center)^2+(all_cells$y-center)^2) > pith_size,]
    all_cells <- rbind(new_cells, all_cells)
    all_cells$id_cell <- 1:nrow(all_cells)

  }

  return(all_cells)
}
