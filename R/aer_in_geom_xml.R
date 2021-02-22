

aer_in_geom_xml <- function(sim, path = "./MECHA/Projects/GRANAR/in/Maize_Geometry.xml"){
  
  require(xml2)
  if (is.null(path)) {
    warning("No path specified")
  }
  if(is.null(sim$id_aerenchyma)){
    warning("No aerenchyma id specified")
  }else{
    id_aerenchyma <- sim$id_aerenchyma
  }
  xml <- read_xml(path)
  
  aer <- xml_children(xml_find_all(xml, "//aerenchyma_range"))
  
  # newbee <- 'aerenchyma id="0"'
  if(length(sim$id_aerenchyma) == 0){
    message("no aerenchyma")
  }else{
    new_siblings <- paste0('aerenchyma id="',id_aerenchyma,'"')
    xml_add_sibling(aer, new_siblings)
  }

  xml_remove(aer[1])
  
  path <- paste0(c(unlist(str_split(path, ".xml"))[1]),"_aer.xml")
  
  write_xml(xml, path)
  return(TRUE)
}