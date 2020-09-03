#' @title Write a input file 'xml' for GRANAR 
#'
#' This function creates a input file 'xml' for GRANAR from dataframe
#' @param path Path to the future XML file containing the different parameters for GRANAR.
#' @param params Table with the different parameters.
#' @author Adrien Heymans and Guillaume Lobet
#' @export
#' @import xml2
#' @examples
#'
#' root_cross_section <- create_anatomy(parameters = param)
#' write_param_xml(path = "Root_id.xml", params = param)


write_param_xml <- function(path, params){
  
  if(!is.null(params$param)){
    colnames(params) <- c("name", "type", "value")
  }
  xml <- "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
  xml <- paste0(xml, "<granar>\n")
  
  
  for(i in unique(params$name)){
    tmp_type <- params$type[params$name == i]
    tmp_val <- params$value[params$name == i]
    xml <- paste0(xml, paste0("\t<",i))
    for (j in 1:length(tmp_type)){
      xml <- paste0(xml," ",tmp_type[j],"=\"",tmp_val[j], "\"")
    }
    xml <- paste0(xml,"\ />\n")
  }
  xml <- paste0(xml, "</granar>")
  if (!is.null(path)) {
    cat(xml, file = path)
    return(TRUE)
  }  else {
    return(xml)
  }
}
