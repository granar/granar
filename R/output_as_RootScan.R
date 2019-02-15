
#'
#' This function turn the GRANAR output into RootScan data.
#' @param sim The output list() generated with create_anatomy() function. Default = NULL
#' @keywords root anatomy RootScan
#' @export
#' @examples
#'
#' data <- output_as_RootScan(sim = create_anatomy(parameters = params))


output_as_RootScan <- function(sim = NULL){
  
  if(is.null(sim)){
    warning("please enter GRANAR output")
    return(NULL)
    }
  
  # Save output data ---------------------------------------
  output <- sim$output
  
  out_data <- tibble(  
    #Total area of the root cross-section
    m_RXSA = output$value[output$io == "output" & output$name == "all" & output$type == "layer_area"]+
      output$value[output$io == "output" & output$name == "aerenchyma" & output$type == "layer_area"],
    m_RXSA_2 = (((max(sim$nodes$x)-min(sim$nodes$x))/2)^2)*pi,
    #Total xylem area
    m_XVA = output$value[output$io == "output" & output$name == "xylem" & output$type == "layer_area"],
    #Total stele area
    m_TSA_2 = (((max(sim$nodes$x[sim$nodes$type == "pericycle"])-min(sim$nodes$x[sim$nodes$type == "pericycle"]))/2)^2)*pi,
    
    m_CCA = output$value[output$io == "output" & output$name == "cortex" & output$type == "layer_area"]+
      output$value[output$io == "output" & output$name == "endodermis" & output$type == "layer_area"]+
      output$value[output$io == "output" & output$name == "exodermis" & output$type == "layer_area"]+
      output$value[output$io == "output" & output$name == "epidermis" & output$type == "layer_area"],
    
    m_CC = output$value[output$io == "output" & output$name == "cortex" & output$type == "n_cells"]+
      output$value[output$io == "output" & output$name == "endodermis" & output$type == "n_cells"]+
      output$value[output$io == "output" & output$name == "exodermis" & output$type == "n_cells"],
    
    m_AA = output$value[output$io == "output" & output$name == "aerenchyma" & output$type == "layer_area"],
    m_pA = output$value[output$io == "output" & output$name == "aerenchyma" & output$type == "proportion"],
    
    granar_time = output$value[output$io == "output" & output$name == "simulation" & output$type == "time"],
    m_TSA = if(output$value[output$io == "input" & output$name == "planttype"] == 1){
      output$value[output$io == "output" & output$name == "stele" & output$type == "layer_area"]+
      output$value[output$io == "output" & output$name == "phloem" &  output$type == "layer_area"]+
      output$value[output$io == "output" & output$name == "pericycle" & output$type == "layer_area"]+
      output$value[output$io == "output" & output$name == "companion_cell" & output$type == "layer_area"]+
      output$value[output$io == "output" & output$name == "xylem" & output$type == "layer_area"]
    }else{
      output$value[output$io == "output" & output$name == "stele" & output$type == "layer_area"]+
      output$value[output$io == "output" & output$name == "pericycle" & output$type == "layer_area"]+
      output$value[output$io == "output" & output$name == "xylem" & output$type == "layer_area"]
    }
    )%>%
    mutate(m_TCA = m_RXSA-m_TSA,
           m_TCA_2 = m_RXSA_2 - m_TSA_2)
  return(out_data)
}