#' @title Merge multiple polygon cells that are adjacent to each other and can form concave polygons.
#'
#'
#' @param data the dataframe containing the cells to merge
#' @keywords root polygons
#' @export
#' @examples
#' # new_cell <- concavety(data = cells)
#'

# merge multiple polygon cells next to each other that can form concave polygons
concavety <- function(data) {

  merged_polygons = list()
  for (i in unique(data$id_cell)) {
    # Get the subset of points for the current id_cell
    subset_points <- rbind(data[data$id_cell == i, ],data[data$id_cell == i, ][1,])

    if (nrow(subset_points) > 2) {
      # Convert points to a matrix for creating an sf object
      poly_matrix <- matrix(c(subset_points$x, subset_points$y), ncol = 2)
      # Create an sf object for the subset points
      subset_sf <- sf::st_polygon(list(poly_matrix))
      # Combine polygons for the current id_cell
      merged_polygon <- sf::st_union(subset_sf)
      # Add the merged polygon to the list
      merged_polygons[[i]] <- merged_polygon
    }
  }

  # Filter out NULL elements from the list
  filtered_list <- Filter(Negate(is.null), merged_polygons)
  # Combine all merged polygons into a single sf object
  polygons_sf <- st_sfc(filtered_list)

  merged_polygons <- st_union(polygons_sf)
  # Clean inner vertices
  cleaned_polygons <- sf::st_simplify(merged_polygons, preserveTopology = FALSE)
  # Extract coordinates as a data frame
  coord <- as.data.frame(sf::st_coordinates(cleaned_polygons)) %>%
    transmute(x = X,
              y = Y)

  new_cell = tibble(id_cell = data$id_cell[1], x = coord$x, y = coord$y, type = data$type[1],
                    area = sf::st_area(cleaned_polygons),
                    dist = mean(unique(data$dist)),
                    angle = mean(unique(data$angle)),
                    radius = mean(unique(data$radius)),
                    id_layer = mean(unique(data$id_layer)),
                    id_group = data$id_group[1],
                    my = mean(coord$y),
                    mx = mean(coord$x),
                    atan = seq(-pi, pi, length.out = length(coord$y)),
                    new = data$id_cell[1],
                    y1 = coord$y,
                    y2 = coord$y[c(2:length(coord$y), 1)],
                    x1 = coord$x,
                    x2 = coord$x[c(2:length(coord$x), 1)],
                    id_point = paste0(coord$x,";",coord$y),
                    euc=NA,
                    out=NA,
                    rank =NA
  )
  return(new_cell%>%select(colnames(data)))
}




