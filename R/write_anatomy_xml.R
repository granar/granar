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
  
  cellgroups <- data.frame(id_group = c(1, 2, 3, 4, 5, 13, 16, 12, 11),
                           type = c("exodermis", "epidermis", "endodermis", "cortex", "stele", "xylem", "pericycle", "companion_cell", "phloem"))
  
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
    
  sim$nodes <- merge(sim$nodes, cellgroups, by="type")  %>% 
    mutate(id_group = id_group.y)
  
  temp_wall <- ddply(sim$nodes, .(id_cell, id_group), summarise, walls = paste0('\t\t\t\t<wall id="', 
                                                                                paste(id_wall-1, collapse='"/>\n\t\t\t\t<wall id="'),
                                                                                '"/>\n'))
  xml <- paste0(xml, paste0('\t\t<cell id="',temp_wall$id_cell-1, '" group="', temp_wall$id_group, '" truncated="false" >\n',
                            '\t\t\t<walls>\n', temp_wall$walls, '\t\t\t</walls>\n',
                            '\t\t</cell>\n', collapse=""))
  xml <- paste0(xml, '\t</cells>\n')
  
  
  # Write the walls information
  xml <- paste0(xml, '\t<walls count="',nrow(sim$walls),'">\n')
  xml <- paste0(xml,paste0('\t\t<wall id="',sim$walls$id_wall-1,'" group="0" edgewall="false" >\n',
                           '\t\t\t<points>\n',
                           '\t\t\t\t<point x="',sim$walls$x1,'" y="',sim$walls$y1,'"/>\n',
                           '\t\t\t\t<point x="',sim$walls$x2,'" y="',sim$walls$y2,'"/>\n',
                           '\t\t\t</points>\n',
                           '\t\t</wall>\n', collapse = ""))
  xml <- paste0(xml, '\t</walls>\n')
  
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