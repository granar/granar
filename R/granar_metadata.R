

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



batch_Gmeta <- function(fls){

  sampl_id = parse_number(fls)

  Out <- tibble(sampl_id = sampl_id,
                RXA = NA, TSA = NA, MXA = NA, AA = NA, TCA = NA,
                TSA_in = NA, MXA_in = NA, nX = NA, kx_unM = NA, kx_M = NA)
  for(j in fls){
    tmp_output <- granar_metadata(paste0("./MECHA/cellsetdata/",j))

    id = parse_number(unlist(str_split(j,"_"))[2])

    TSA_in <- pi*(tmp_output$value[tmp_output$io == "input" & tmp_output$param == "stele" & tmp_output$type == "layer_diameter"]/2)^2

    MXA_in <- tmp_output$value[tmp_output$io == "input" & tmp_output$param == "xylem" & tmp_output$type == "n_files"]*pi*(tmp_output$value[tmp_output$io == "input" & tmp_output$param == "stele" & tmp_output$type == "cell_diameter"]+ tmp_output$value[tmp_output$io == "input" & tmp_output$param == "xylem" & tmp_output$type == "max_size"]/2)^2
    TSA <- tmp_output$value[tmp_output$io == "output" & tmp_output$param == "stele" & tmp_output$type == "layer_area"]+tmp_output$value[tmp_output$io == "output" & tmp_output$param == "pericycle" & tmp_output$type == "layer_area"]
    MXA <- tmp_output$value[tmp_output$io == "output" & tmp_output$param == "xylem" & tmp_output$type == "layer_area"]
    RXA <- tmp_output$value[tmp_output$io == "output" & tmp_output$param == "all" & tmp_output$type == "layer_area"]
    TCA <- tmp_output$value[tmp_output$io == "output" & tmp_output$param == "cortex_to_epidermis" & tmp_output$type == "layer_area"]
    AA <- tmp_output$value[tmp_output$io == "output" & tmp_output$param == "aerenchyma" & tmp_output$type == "layer_area"]
    if(length(AA) == 0){AA = 0}

    nX <- tmp_output$value[tmp_output$io == "input" & tmp_output$param == "xylem" & tmp_output$type == "n_files"]
    PXA_1 <- tmp_output$value[tmp_output$io == "output" & tmp_output$param == "stele" & tmp_output$type == "cell_area"]
    ratio <- tmp_output$value[tmp_output$io == "input" & tmp_output$param == "xylem" & tmp_output$type == "ratio"]
    nPX = nX*ratio
    PXA = nPX*PXA_1
    a = PXA_1*1000^2
    k_protxyl_s = a^2/(8*pi*200*1E-5/3600/24)*1E-12
    # kx when only the proto xylem have their cell wall lignified
    kx_unM = k_protxyl_s*nPX*200/1E4

    LMXA = MXA - PXA
    LMXA_1 = LMXA/nX
    a = LMXA_1*1000^2
    k_Mxyl_s = a^2/(8*pi*200*1E-5/3600/24)*1E-12
    # kx when all xylem elements have their cell wall lignified
    kx_M = k_Mxyl_s*nX*200/1E4 + kx_unM

    Out$RXA[Out$sampl_id == id] = RXA
    Out$TSA[Out$sampl_id == id] = TSA
    Out$MXA[Out$sampl_id == id] = MXA
    Out$TCA[Out$sampl_id == id] = TCA
    Out$AA[Out$sampl_id == id] = AA
    Out$TSA_in[Out$sampl_id == id] = TSA_in
    Out$MXA_in[Out$sampl_id == id] = MXA_in
    Out$nX[Out$sampl_id == id] = nX
    Out$kx_M[Out$sampl_id == id] = kx_M
    Out$kx_unM[Out$sampl_id == id] = kx_unM
  }

  return(Out)
}

