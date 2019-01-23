

#' Get GRANAR metadata out of .xml
#'
#' To get the metadata out of the xml file created with write_anatomy_xml.
#' @param path Path to the XML file containing granar output. Default = NULL.
#' @keywords root anatomy
#' @export
#' @examples
#' metadata <- granar_metadata("PATH TO THE XML FILE")
#'


granar_metadata <- function(path = NULL){
  
  if(is.null(path) ){warning("No path specified")}
  x <- read_xml(path)
  
  # GET THE CELLS DATA
  mydata <- xml_children(x)[2]
  temp <- xml_find_all(x, ".//metadata")
  param <- xml_find_all(temp, ".//parameter")
  output <- NULL
  for(i in c(1:length(param))){
    output <- rbind(output, data.frame(io = xml_attr(param[i], "io"),
                                     param = xml_attr(param[i], "name"),
                                     type = xml_attr(param[i], "type"),
                                     value = as.numeric(xml_attr(param[i], "value"))))
  }
  return(output)
}
  