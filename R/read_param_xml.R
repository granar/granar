#' Read the parameters for GRANAR
#'
#' Read the parameters for GRANAR from an XML file
#' @param path The path to the XML files ith the parameters
#' @keywords granar
#' @export
#' @import xml2
#' @examples
#' read_param_xml()
#'

read_param_xml <- function(path = NULL){

  if( is.null(path) ){
    warning("No path specified")
  }

  input <- read_xml(path)
  params <- NULL

  # Quality checks. Check if all the needed tags are present in the XML file
  to_find <- c("planttype", "randomness", "xylem", "phloem", "stele", "endodermis", "exodermis", "epidermis", "aerenchyma", "pericycle", "cortex")
  for(tf in to_find){
    if (length(xml_find_all(input, paste0("//",tf))) == 0) warning(paste0("Could not find the '",tf,"' tag in the XML file"))
  }

  #Read the file and get the parameters in a table
  for( ch in xml_children(xml_find_all(input, "//*"))){
    att <- xml_attrs(ch)
    for(i in c(1:length(att))){
      params <- rbind(params, data.frame(
        name = xml_name(ch),
        type = names(att)[i],
        value = att[i]
      ))
    }
  }
  row.names(params) <- NULL
  params <- params %>% mutate(value = as.numeric(as.character(value)))

  return(params)
}
