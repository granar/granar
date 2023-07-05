#' @title Prepare cell wall data
#'
#' tydy cells as vertex data
#' @param rs1 The nodes dataframe
#' @keywords root
#' @export
#'

vertex <- function(rs1){
  nodes <- rs1 %>%
    mutate(id_point= paste0(rs1$x,";",rs1$y))%>%
    group_by(id_cell) %>%
    filter(!duplicated(id_point))%>%
    dplyr::mutate(xx = c(x[-1],x[1])) %>%
    dplyr::mutate(yy = c(y[-1],y[1]))%>%
    select(-id_point)

  nodes <- nodes %>%
    ungroup() %>%
    mutate(vertical = ifelse(x == xx, "true", "false")) %>%
    mutate(x1 = ifelse(x > xx, x, xx)) %>%
    mutate(x2 = ifelse(x > xx, xx, x)) %>%
    mutate(y1 = ifelse(x > xx, y, yy)) %>%
    mutate(y2 = ifelse(x > xx, yy, y)) %>%
    # If wall is perfectly vertical
    mutate(y1 = ifelse(x == xx,
                       ifelse(y > yy, yy, y), y1)) %>%
    mutate(y2 = ifelse(x == xx,
                       ifelse(y > yy, y, yy), y2)) %>%
    mutate(wall_length = sqrt((x2-x1)^2 + (y2-y1)^2)) %>%
    mutate(wall_length2 = sqrt((xx-x)^2 + (yy-y)^2),
           slope = (y2-y1)/(x2-x1),
           intercept = y1 - slope*x1)
  return(nodes)
}
