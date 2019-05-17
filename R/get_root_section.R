#' Read GRANAR file under XML
#'
#' This function load the GRANAR file .XML.
#' @param path The path to the GRANAR file .XML
#' @keywords root
#' @export
#' @examples
#' root <- list()
#' root$nodes <- get_root_section("PATH_TO_XLM_FILE")
#' plot_anatomy(root, col = "type")

get_root_section <- function(path){

  if(is.null(path) ){warning("No path specified")}
  x <- read_xml(path)

  # GET THE CELLS DATA
  mydata <- xml_children(x)[2]
  temp <- xml_find_all(x, ".//cell")
  cells <- NULL
  for(i in c(1:length(temp))){
    wall <- xml_find_all(temp[i], ".//wall")
    n <- length(wall)
    cells <- rbind(cells, data.frame(id_cell = rep(as.numeric(xml_attr(temp[i], "id")), n),
                                     group = rep(as.numeric(xml_attr(temp[i], "group")), n),
                                     wall_id = as.numeric(xml_attr(wall, "id"))))
  }

  # GET THE WALL DATA
  mydata <- xml_children(x)[3]
  temp <- xml_find_all(mydata, ".//wall")
  walls <- NULL
  for(i in c(1:length(temp))){
    points <- xml_find_all(temp[i], ".//point")
    n <- length(points)
    walls <- rbind(walls, data.frame(id_cell = rep(as.numeric(xml_attr(temp[i], "id")), n),
                                     group = rep(as.numeric(xml_attr(temp[i], "group")), n),
                                     x = as.numeric(xml_attr(points, "x")),
                                     y = as.numeric(xml_attr(points, "y"))))
  }


  # GET THE GROUP INFORMATIONS
  mydata <- xml_find_all(x, ".//group")
  groups <- data.frame(id_cell = as.numeric(xml_attr(mydata, "id")),
                       type = xml_attr(mydata, "name")
  )

  #MERGE THE DATA
  rs <- merge(cells, walls, by.x = "wall_id", by.y = "id_cell")
  rs <- merge(rs, groups, by.x = "group.x", by.y = "id_cell")

  # REORDER THE POINTS WITHIN EACH CELL tO HAVE A NICE POLYGON
  rs1 <- NULL
  for(i in unique(rs$id_cell)){
    temp <- rs[rs$id_cell == i,]
    temp$atan <- atan2(temp$y - mean(temp$y), temp$x - mean(temp$x))
    temp <- temp[order(temp$atan),]
    temp <- temp[!duplicated(temp$atan),]
    rs1 <- rbind(rs1,temp)
  }
  root_section <- rs1[,c("id_cell", "x", "y", "type")]
  return(root_section)
}

