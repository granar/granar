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
  
  cellgroups <- data.frame(id_group = c(1, 2, 3, 4, 5, 13, 16),
                           type = c("exodermis", "epidermis", "endodermis", "cortex", "stele", "xylem", "pericycle"))
  
  xml <- '<?xml version="1.0" encoding="utf-8"?>\n'
  xml <- paste0(xml, '<crosssimdata>\n')
  
  # Metadata
  xml <- paste0(xml, '\t<metadata>\n')
  xml <- paste0(xml, '\t\t<parameters>\n')
  xml <- paste0(xml, '\t\t\t<parameter name="num_cortex" value="',num_cortex,'"/>\n')
  xml <- paste0(xml, '\t\t\t<parameter name="diam_cortex" value="',diam_cortex,'"/>\n')
  xml <- paste0(xml, '\t\t\t<parameter name="size_stele" value="',size_stele,'"/>\n')
  xml <- paste0(xml, '\t\t\t<parameter name="diam_stele" value="',diam_stele,'"/>\n')
  xml <- paste0(xml, '\t\t\t<parameter name="proportion_aerenchyma" value="',proportion_aerenchyma/100,'"/>\n')
  xml <- paste0(xml, '\t\t\t<parameter name="n_aerenchyma_files" value="',n_aerenchyma_files,'"/>\n')
  xml <- paste0(xml, '\t\t\t<parameter name="n_xylem_files" value="',n_xylem_files,'"/>\n')
  xml <- paste0(xml, '\t\t\t<parameter name="diam_xylem" value="',diam_xylem,'"/>\n')
  xml <- paste0(xml, '\t\t</parameters>\n')
  xml <- paste0(xml, '\t</metadata>\n')
  
  
  # Cells
  xml <- paste0(xml, '\t<cells count="',nrow(sim$cells),'">\n')
  sim$nodes <-merge(sim$nodes, cellgroups, by="type")
  temp_wall <- ddply(sim$nodes, .(id_cell, id_group), summarise, walls = paste0('\t\t\t\t<wall id="', 
                                                                                paste(id_wall-1, collapse='"/>\n\t\t\t\t<wall id="'),
                                                                                '"/>\n'))
  xml <- paste0(xml, paste0('\t\t<cell id="',temp_wall$id_cell-1, '" group="', temp_wall$id_group, '" truncated="false" >\n',
                            '\t\t\t<walls>\n', temp_wall$walls, '\t\t\t</walls>\n',
                            '\t\t</cell>\n', collapse=""))
  xml <- paste0(xml, '\t</cells>\n')
  
  
  # Walls
  xml <- paste0(xml, '\t<walls count="',nrow(sim$walls),'">\n')
  xml <- paste0(xml,paste0('\t\t<wall id="',sim$walls$id_wall-1,'" group="0" edgewall="false" >\n',
                           '\t\t\t<points>\n',
                           '\t\t\t\t<point x="',sim$walls$x1,'" y="',sim$walls$y1,'"/>\n',
                           '\t\t\t\t<point x="',sim$walls$x2,'" y="',sim$walls$y2,'"/>\n',
                           '\t\t\t</points>\n',
                           '\t\t</wall>\n', collapse = ""))
  xml <- paste0(xml, '\t</walls>\n')
  
  # Groups
  xml <- paste0(xml, '\t<groups>\n')
  xml <- paste0(xml, '\t\t<cellgroups>\n')
  for(i in c(1:nrow(cellgroups))){
    xml <- paste0(xml, '\t\t\t<group id="',cellgroups$id[i],'" name="',cellgroups$name[i],'" />\n')
  }
  xml <- paste0(xml, '\t\t</cellgroups>\n')
  xml <- paste0(xml, '\t\t<wallgroups>\n')
  xml <- paste0(xml, '\t\t\t<group id="0" name="unassigned" />\n')
  xml <- paste0(xml, '\t\t</wallgroups>\n')
  xml <- paste0(xml, '\t</groups>\n')
  
  xml <- paste0(xml, '</crosssimdata>')
  
  if(!is.null(path)){
    cat(xml, file = path)
    return(TRUE)
  }else{
    return(xml)
  }
  
  
}