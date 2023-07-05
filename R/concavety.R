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

    # Create SpatialPointsDataFrame from x and y coordinates
    points <- sp::SpatialPointsDataFrame(
      data.frame(x = data$x, y = data$y),
      data = data,
      proj4string = CRS("+proj=longlat +datum=WGS84")
    )

    # Create a list to store the polygons
    merged_polygons <- list()

    for (i in unique(points$id_cell)) {
      # Get the subset of points for the current id_cell
      subset_points <- points[points$id_cell == i, ]

      if(nrow(subset_points)> 2){
        # Convert points to polygons
        polygons <- suppressWarnings(SpatialPolygons(list(
          Polygons(list(Polygon(coordinates(subset_points))), ID = as.character(i)))))
        # Merge the polygons for the current id_cell
        merged_polygon <- maptools::unionSpatialPolygons(polygons, IDs = i)
        # Add the merged polygon to the list
        merged_polygons[[i]] <- merged_polygon
      }
    }
    filtered_list <- Filter(Negate(is.null), merged_polygons)
    # Combine all merged polygons into a single SpatialPolygonsDataFrame
    merged_polygons <- do.call(rbind, filtered_list)
    merged_polygons_df <- SpatialPolygonsDataFrame(
      merged_polygons,
      data = data.frame(id_cell = unique(data$id_cell), id_group = data$id_group[1]),
      match.ID = FALSE
    )
    # Remove inner veritiles
    merged_sp_clean <- rgeos::gUnaryUnion(merged_polygons_df)
    poly_clean = merged_sp_clean@polygons[[1]]
    coord = as.data.frame( poly_clean@Polygons[[1]]@coords)
    new_cell = tibble(id_cell = data$id_cell[1], x = coord$x, y = coord$y, type = data$type[1],
                      area = poly_clean@Polygons[[1]]@area,
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

