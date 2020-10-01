#' Write the root anatomy as an XML file.
#'
#' The structure of the XL file matches the one of CellSet
#' @param sim The simulation results
#' @param path The path where to save the xml file
#' @keywords root
#' @export
#' @examples
#' write_anatomy_xml()
#'

write_anatomy_xml <- function(sim = NULL, path = NULL){

  if(is.null(sim)) warning("No simulation found. Please input a GRANAR simulation")
  if(is.null(path)) warning("No path found to save the XML file")


  if(length(sim$walls$x3) > 0){
    nodal <- sim$walls_nodes
  }else{
    nodal <- sim$nodes
  }
  nodal <- nodal%>%
    select(-id_group)

  cellgroups <- data.frame(id_group = c(1, 2, 3, 3, 4, 5, 13, 16, 12, 11, 4, 4, 3),
                           type = c("exodermis", "epidermis", "endodermis", "passage_cell",  "cortex", "stele", "xylem", "pericycle", "companion_cell", "phloem", "inter_cellular_space", "aerenchyma", "cambium"))

  xml <- '<?xml version="1.0" encoding="utf-8"?>\n'
  xml <- paste0(xml, '<granardata>\n')

  # Write the Metadata
  xml <- paste0(xml, '\t<metadata>\n')
  xml <- paste0(xml, '\t\t<parameters>\n')
  xml <- paste0(xml,paste0('\t\t\t<parameter io="',sim$output$io,'" ',
                            'name="',sim$output$name,'" ',
                            'type="',sim$output$type,'" ',
                            'value="',sim$output$value,'"/>\n', collapse = ""))
  xml <- paste0(xml, '\t\t</parameters>\n')
  xml <- paste0(xml, '\t</metadata>\n')

  # Write the cells information
  xml <- paste0(xml, '\t<cells count="',nrow(sim$cells),'">\n')

  nodes_data <- merge(nodal, cellgroups, by="type")

  temp_wall <- ddply(nodes_data, .(id_cell, id_group), summarise, walls = paste0('\t\t\t\t<wall id="',
                                                                                paste(id_wall-1, collapse='"/>\n\t\t\t\t<wall id="'),
                                                                                '"/>\n'))
  xml <- paste0(xml, paste0('\t\t<cell id="',temp_wall$id_cell-1, '" group="', temp_wall$id_group, '" truncated="false" >\n',
                            '\t\t\t<walls>\n', temp_wall$walls, '\t\t\t</walls>\n',
                            '\t\t</cell>\n', collapse=""))
  xml <- paste0(xml, '\t</cells>\n')


  # Write the walls information
  xml <- paste0(xml, '\t<walls count="',length(unique(nodes_data$id_wall)),'">\n')

  walls <- sim$walls%>%
    select((starts_with("x") | starts_with("y")) & ends_with(as.character(c(0:9))))
  col_nam <- colnames(walls)

  substr1(col_nam[nchar(col_nam) == 2], 1.5) <- "0"
  col_nam <- paste0(substr1(col_nam, -1), "_", substr1(col_nam, 1))
  colnames(walls) <- col_nam

  sorted_name <- sort(col_nam)
  walls <- walls%>%select(sorted_name)

  N <- max(parse_number(col_nam))
  begin <- tibble(tag1 = '\t\t<wall id="',
                id_wall = sim$walls$id_wall-1,
                tag2 = '" group="0" edgewall="false" >\n\t\t\t<points>\n')
  middle <- tibble(tag_x1 = '\t\t\t\t<point x="',
                   x1 = walls[,sorted_name[1]],
                   tag_y1 = '" y="',
                   y1 = walls[,sorted_name[2]],
                   tag_end1 = '"/>\n')
  for(k in 2:N){
    h <- k*2-1 # odd number
    tmp_coord <- walls[,sorted_name[c(h,h+1)]]
    tmp_middle <- tibble(tag_x = '\t\t\t\t<point x="',
                         x = tmp_coord[,1],
                         tag_y = '" y="',
                         y = tmp_coord[,2],
                         tag_end = '"/>\n')
    tmp_col_name <- colnames(tmp_middle)
    colnames(tmp_middle) <- paste0(t(tmp_col_name), k)
    middle <- cbind(middle, tmp_middle)
  }
  taged_walls <- cbind(begin,middle)%>%
    mutate(tag_ending = '\t\t\t</points>\n\t\t</wall>\n')
  xml <- paste0(xml, paste0(t(taged_walls), collapse = ""))
  xml <- paste0(xml, '\t</walls>\n')
  xml <- str_remove_all(xml, '\t\t\t\t<point x=\"NA\" y=\"NA\"/>\n')

  # Write the cell group informations
  print(cellgroups)
  xml <- paste0(xml, '\t<groups>\n')
  xml <- paste0(xml, '\t\t<cellgroups>\n')
  for(i in c(1:nrow(cellgroups))){
    xml <- paste0(xml, '\t\t\t<group id="',cellgroups$id_group[i],'" name="',cellgroups$type[i],'" />\n')
  }
  xml <- paste0(xml, '\t\t</cellgroups>\n')
  xml <- paste0(xml, '\t\t<wallgroups>\n')
  xml <- paste0(xml, '\t\t\t<group id="0" name="unassigned" />\n')
  xml <- paste0(xml, '\t\t</wallgroups>\n')
  xml <- paste0(xml, '\t</groups>\n')

  xml <- paste0(xml, '</granardata>')

  if(!is.null(path)){
    cat(xml, file = path)
    return(TRUE)
  }else{
    return(xml)
  }


}


substr1 <- function(x,y) {
  z <- sapply(strsplit(as.character(x),''),function(w) paste(na.omit(w[y]),collapse=''))
  dim(z) <- dim(x)
  return(z) }

`substr1<-` <- function(x,y,value) {
  names(y) <- c(value,rep('',length(y)-length(value)))
  z <- sapply(strsplit(as.character(x),''),function(w) {
    v <- seq(w)
    names(v) <- w
    paste(names(sort(c(y,v[setdiff(v,y)]))),collapse='') })
  dim(z) <- dim(x)
  return(z) }
