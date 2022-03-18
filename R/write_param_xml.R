#' @title Write the parameters of the root anatomy as an XML file.
#'
#' The structure of the XL file matches the default input files of GRANAR
#' @param sim The simulation results
#' @param path The path where to save the xml file
#' @keywords root
#' @export
#' @examples
#' write_param_xml(sim, path = "params.xml" )
#'

write_param_xml <- function(sim = NULL, path = NULL){

  if(is.null(sim)) warning("No simulation found. Please input a GRANAR simulation")
  if(is.null(path)) warning("No path found to save the XML file")

  param <- sim$output%>%
      filter(io == "input")

  xml <- '<?xml version="1.0" encoding="utf-8"?>\n'
  xml <- paste0(xml, '<granar>\n')


  for(i in unique(param$name)){

    xml <- paste0(xml, '<',i)
    for(j in param$type[param$name == i]){
      xml <- paste0(xml, ' ',j, '="',param$value[param$name == i & param$type == j],'"')
    }
    xml <- paste0(xml, '/>\n')
  }

  xml <- paste0(xml, '</granar>')

  if(!is.null(path)){
    cat(xml, file = path)
    return(TRUE)
  }else{
    return(xml)
  }

}
