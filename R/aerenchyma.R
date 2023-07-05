
#' @title Create aerenchyma
#'
#'
#' @param params input parameter
#' @param rs1 node dataframe
#' @keywords root
#' @export
#' @examples
#' # rs <- aerenchyma(params, rs1)
#'

aerenchyma <- function(params, rs1){

  n_aerenchyma_files <- params$value[params$name == "aerenchyma" & params$type == "n_files"]
  proportion_aerenchyma <- params$value[params$name == "aerenchyma" & params$type == "proportion"]

  angle_inc <- (2 * pi) / n_aerenchyma_files
  safe_cortex_layer <- c(min(rs1$id_layer[rs1$type == "cortex"]))
  last_cortex_layer <- c(max(rs1$id_layer[rs1$type == "cortex"]))
  area_all <- rs1%>%
    filter(type %in% c("cortex", "inter_cellular_space", "exodermis", "epidermis"))%>%
    dplyr::group_by(id_cell)%>%
    dplyr::summarise(area = mean(area))
  ini_cortex_area <- sum(area_all$area)
  surface_to_kill <- ini_cortex_area*proportion_aerenchyma
  stk_zone = surface_to_kill/n_aerenchyma_files


  # select type of aerenchyma formation
  aer_type <- params$value[params$name == "aerenchyma" & params$type == "type"]
  if(length(aer_type)== 0){
    aer_type <- params$value[params$name == "planttype" & params$type == "param"]
  }
  if (aer_type == 1){
    small_r <- mean(rs1$dist[rs1$id_layer == safe_cortex_layer])+0.5*mean(rs1$radius[rs1$id_layer == safe_cortex_layer])
    big_R <- mean(rs1$dist[rs1$id_layer == last_cortex_layer])+0.5*mean(rs1$radius[rs1$id_layer == last_cortex_layer])
    angle_range_inc <- stk_zone/(big_R^2-small_r^2)

  }else if(aer_type == 2){
    angle_range_inc <- (2 * pi * proportion_aerenchyma / 100) / n_aerenchyma_files
  }else{warning("aerenchyma type is out of range [1-2]" )}

  angle_range_inc_ini <- angle_range_inc
  cortex_area <- ini_cortex_area
  angle <- stats::runif(1, 0.6, 1) * pi/n_aerenchyma_files
  angle_range <- c(angle - angle_range_inc, angle + angle_range_inc)
  id_group_max <- max(rs1$id_group)

  for(j in c(1:n_aerenchyma_files)){
    # saved_rs1 <- rs1
    possi <- rs1[rs1$id_layer %!in% safe_cortex_layer & rs1$type %in% c("cortex", "inter_cellular_space"),]
    gotogo <- T
    n_try <- 0
    while(gotogo){
      angle_range <- c(angle - angle_range_inc, angle + angle_range_inc)
      ma_cell <- rs1 %>%
        filter(id_layer %!in% safe_cortex_layer ,
               type == "cortex" ,
               angle > angle_range[1] ,
               angle < angle_range[2])
      are <- sum_area(ma_cell)
       if(are > stk_zone){
        angle_range_inc <- angle_range_inc*0.9
        n_try <- n_try + 1
        if(n_try > 100){
          gotogo <- F
        }
        next()
      }

      if(nrow(ma_cell) > 0){
        while(are <= stk_zone ){

          are <- sum_area(ma_cell)

          pointy <- unique(ma_cell$id_point)
          done <- unique(ma_cell$id_cell) # id from potential cell to be named aer
          nei <- unique(possi$id_cell[possi$id_cell %!in% done & possi$id_point %in% pointy]) # neigh cells and inter_cell
          icp <- unique(possi$id_cell[possi$id_cell %in% nei & possi$type == "inter_cellular_space"]) # neigh inter_cell
          pre <- unique(possi$id_cell[possi$id_cell %in% c(done,icp)]) # merge potential and inter_cell

          done <- c(done, nei) # full potential
          ma_cell <- possi[possi$id_cell %in% done,]

          are <- sum_area(ma_cell)

          n_try <- n_try + 1
          if(n_try > 100){
            are = stk_zone
          }

        }
        # once it overflow the threshold limit

        ma_cell <- possi[possi$id_cell %in% pre,] # this default potential cell and inter_cell are selected

        are_b <- sum_area(ma_cell)

        # adjust to include the right amount of cells as aerenchyma
        if(are_b > stk_zone){to_be_added <- NULL}else{

          next_area <- possi[possi$id_cell %in% done & possi$id_cell %!in% pre,]

          alm <- next_area%>%
            dplyr::group_by(id_cell)%>%
            dplyr::summarise(area = mean(area))%>%
            ungroup()

          alm = alm%>%
            mutate(sarea = cumsum(area),
                   near = abs(sarea+are_b-stk_zone))

          jk <- which(alm$near == min(alm$near))

          to_be_added <- alm$id_cell[1:jk]
          ma_cell <- possi[possi$id_cell %in% c(pre,to_be_added),]
          pointy <- unique(ma_cell$id_point)
          done <- unique(ma_cell$id_cell) # id from potential cell to be named aer
          nei <- unique(possi$id_cell[possi$id_cell %!in% done & possi$id_point %in% pointy]) # neigh cells and inter_cell
          icp <- unique(possi$id_cell[possi$id_cell %in% nei & possi$type == "inter_cellular_space"])
          pre <- unique(possi$id_cell[possi$id_cell %in% c(done,icp)])
          ma_cell <- possi[possi$id_cell %in% pre,]
        }

        rs1 <- rs1[rs1$id_cell %!in% pre,]
        gotogo <- F # we are done for this part
        #ma_cell = saved_ma_cell <- ma_cell
        aer <- concavety(ma_cell)
        aer$type <- as.character("aerenchyma")
        aer$id_group <- j+id_group_max


        rs1 <- rbind(rs1, aer)

      }else{
        angle_range_inc <- angle_range_inc*1.1
        n_try <- n_try + 1
        if(n_try > 100){
          gotogo <- F
        }
      }

    }
    angle <- angle + angle_inc
    message(paste0(j,"/",n_aerenchyma_files))
    # j <- j +1
  }
  rs <- rs1

  return(rs)

}
