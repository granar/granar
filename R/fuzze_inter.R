#' @title Merge connected inter-cellular space together
#'
#' @param rs1 node data frame
#' @keywords root
#' @export
#'

fuzze_inter <- function(rs1){
  # saved_rs <- rs1
  space <- rs1[rs1$type == "inter_cellular_space",]
  if(nrow(space)> 0){

    double <- space%>%
      dplyr::group_by(id_point)%>%
      dplyr::summarise(n = n())%>%
      filter(n > 1)
    to_correct <- unique(space$id_cell[space$id_point %in% double$id_point])

    if(nrow(double) > 0){
      done <- NULL
      itm <- 0
      for (btw in to_correct) {
        comu1 <- space$id_point[space$id_cell == btw & space$id_point %in% double$id_point]
        nei <- unique(space$id_cell[space$id_cell != btw & space$id_point %in% comu1])
        bou <- c(btw, nei)
        ke <- length(bou)
        te <- 0
        while(ke > te){
          te <- length(bou)
          comu1 <- space$id_point[space$id_cell %in% bou & space$id_point %in% double$id_point]
          nei <- unique(space$id_cell[space$id_cell %!in% bou & space$id_point %in% comu1])
          bou <- unique(c(bou, nei))
          ke <- length(bou)
        }
        itm <- bou

        if(itm[1] %in% done){next()}
        tmp_cell <- space[space$id_cell %in% itm, ]
        if(nrow(tmp_cell) > 0){
          tmp_cell <- tmp_cell%>%
            dplyr::group_by(id_cell)%>%
            dplyr::filter(!duplicated(id_point))

          rs1 <- rs1[rs1$id_cell %!in% itm,]
          tmp_cell <- concavety(tmp_cell)
          rs1 <- rbind(rs1,tmp_cell)
        }else{
          rs1 <- rs1[rs1$id_cell %!in% itm,]
        }
        done <- c(done, itm)
      }

    }
  }

  return(rs1)
}
