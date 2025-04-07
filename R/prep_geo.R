

#' @title Pre-proc for .geo file generation
#'
#' Add cell wall thickness
#' smooth cell corners
#' geo file can be use in GMSH (Require GMSH https://gmsh.info/)
#' @param cross_section The nodes dataframe of the cross section, the coordinates should be in micron
#' @param cell_wall_thickness the inner wall thickness data frame per cell type
#' @param corner_smoothing smoothing algorithm for the cells, can be cell type specific
#' @keywords cell wall GMSH Geo
#' @import smoothr
#' @import sf
#' @export
#'
#'

prep_geo <- function(cross_section, cell_wall_thickness = 0.2, corner_smoothing=0.5){

  if(length(cell_wall_thickness)==1){
    cell_wall_thickness = tibble(type = "default", value = cell_wall_thickness)
  }
  if(length(cell_wall_thickness$value[cell_wall_thickness$type == "outerwall"])==0){
    cell_wall_thickness = rbind(cell_wall_thickness,
                                tibble(type = "outerwall",
                                       value = cell_wall_thickness$value[cell_wall_thickness$type == "default"]))
  }

  if(length(corner_smoothing)==1){
    corner_smoothing = tibble(type = "default", value = corner_smoothing)
  }

  id_cell_vector = unique(cross_section$id_cell)
  root_cell =NULL
  polygons_list = list()
  k = 1
  for(i in id_cell_vector){

    tmp = cross_section%>% filter(id_cell == i)
    wall_type = tmp$type[1]
    if(wall_type %!in% cell_wall_thickness$type){wall_type = "default"}

    sf_linestring <- sf::st_sfc(st_linestring(as.matrix(rbind(tmp[, c("x", "y")],tmp[1, c("x", "y")]))), crs = 2056)
    my_multilinestring = sf::st_sf(geom = sf::st_sfc(sf_linestring), crs = 2056)
    r_poly <-  sf::st_union(my_multilinestring)%>% sf::st_polygonize() %>% sf::st_collection_extract()
    r_poly_smooth <- smoothr::smooth(r_poly, method = "ksmooth",
                                     smoothness = corner_smoothing$value[corner_smoothing$type == wall_type][1])
    shrunken_polygon <- st_buffer(r_poly_smooth,
                                  -cell_wall_thickness$value[cell_wall_thickness$type == wall_type])
    shrunken_polygon  <- st_simplify(shrunken_polygon, dTolerance = 0.1, preserveTopology = TRUE)

    swollen_polygon <- st_buffer(r_poly, cell_wall_thickness$value[cell_wall_thickness$type == "outerwall"])
    # Get the coordinates
    coords <- sf::st_coordinates(shrunken_polygon )[, 1:2]

    pol = tibble(x = coords[,1], y = coords[,2],
                 id_cell = i)
    pol = pol %>%
      mutate(x1 = x,
             y1 = y,
             x2 = c(pol$x[-1], pol$x[1]),
             y2 = c(pol$y[-1], pol$y[1]))

    root_cell = rbind(root_cell, pol)

    polygons_list[[k]] = swollen_polygon
    k = k+1
  }

  # Convert list to sf object
  polygons_sf <- do.call(rbind, lapply(seq_along(polygons_list), function(i) {
    st_sf(id = i, geometry = st_geometry(polygons_list[[i]]))
  }))

  polygons_sfc <- st_as_sfc(polygons_sf)
  # Perform the union operation
  final_polygon <- st_union(polygons_sfc)

  if(length(corner_smoothing$value[corner_smoothing$type == "outerwall"])==0){
    corner_smoothing = rbind(corner_smoothing, tibble(type = "outerwall",
                                                      value = corner_smoothing$value[corner_smoothing$type == "default"]))
  }
  final_polygon <- smoothr::smooth(final_polygon, method = "ksmooth",
                                   smoothness = corner_smoothing$value[corner_smoothing$type == "outerwall"])
  simplified_polygon <- st_simplify(final_polygon, dTolerance = 0.1, preserveTopology = TRUE)
  # Get the coordinates
  wall_coords <- sf::st_coordinates(simplified_polygon)[, 1:2]

  wall = tibble(x = wall_coords[,1], y = wall_coords[,2], id_cell = k)
  wall = wall %>%
    mutate(x1 = x,
           y1 = y,
           x2 = c(wall$x[-1], wall$x[1]),
           y2 = c(wall$y[-1], wall$y[1]))

  return(rbind (wall, root_cell)%>%mutate(res = 1))
}

`%!in%` <- compose(`!`, `%in%`)



