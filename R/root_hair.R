#' @title Create root hair from epidermis cells
#'
#' extend some cells (selected randomly) to make root hair
#' @param params The input dataframe
#' @param rs1 The nodes dataframe
#' @param center The cross-section center
#' @keywords root
#' @export
#'
#'


root_hair <- function(rs1, params, center){

  # initialize variable
  random_fact <- params$value[params$name == "randomness"] / 10 * params$value[params$name == "stele" & params$type == "cell_diameter"]
  n_hair <- params$value[params$name == "hair" & params$type == "n_files"]
  len_hair <- params$value[params$name == "hair" & params$type == "length"]
  hair_r <- params$value[params$name == "epidermis" & params$type == "cell_diameter"]/2

  # existing contition
  if(length(n_hair) != 0){
    x0 <- center
    y0 <- center
    # Select randomly which cell has a root hair
    id_h <- sample(unique(rs1$id_cell[rs1$type == "epidermis"]), n_hair)
    hairy <- NULL

    # For each root hair cell
    for (h in id_h) {

      tmp_cell <- rs1%>%filter(id_cell == h)
      tmp_cell <- tmp_cell%>%
        mutate(euc = sqrt((x0-x)^2 + (y0-y)^2),
               out = ifelse(euc > mean(euc), "out", "in"), # select point in the outer part of the cell
               rank = rank(euc))
      # if the epidermis cell has more than two points (a membrane and not a point) touching the outside
      if(length(tmp_cell$out[tmp_cell$out == "out"])>2){

        junc <- tmp_cell[tmp_cell$out == "out" & tmp_cell$euc != max(tmp_cell$euc),]

        # draw line in between the center of the cross section and the cell mass center
        x1 <- tmp_cell$mx[1]
        y1 <- tmp_cell$my[1]
        m <- (y1-y0)/(x1-x0)
        d <- y0-m*x0
        # Length of the root hair
        l_hair <- len_hair+ len_hair*stats::runif(1,-1,1)*random_fact*100
        sig <- (l_hair^2)*(1+m^2)-(y1-m*x1-d)^2
        sigy <- ((l_hair-hair_r)^2)*(1+m^2)-(y1-m*x1-d)^2
        x2_p <- (x1+y1*m-d*m+sqrt(sig))/(1+m^2)
        x2_m <- (x1+y1*m-d*m-sqrt(sig))/(1+m^2)
        y2_p <- (d+x1*m+y1*m^2+sqrt(sig)*m)/(1+m^2)
        y2_m <- (d+x1*m+y1*m^2-sqrt(sig)*m)/(1+m^2)
        di <- rbind(tibble(euc = sqrt((x0-x2_p)^2+(y0-y2_p)^2), pm = "p"),
                    tibble(euc = sqrt((x0-x2_m)^2+(y0-y2_m)^2), pm = "m"))

        if(di$pm[di$euc == max(di$euc)] == "p"){

          #central point of the root hair apex
          x3 <- (x1+y1*m-d*m+sqrt(sigy))/(1+m^2)
          y3 <- (d+x1*m+y1*m^2+sqrt(sigy)*m)/(1+m^2)
          # perpendicular line
          inv_m <- -1/m
          d3 <- y3-inv_m*x3
          # distance
          s <- (hair_r^2)*(1+inv_m^2)-(y3-inv_m*x3-d3)^2
          #coordinate of the point
          x4 <- (x3+y3*inv_m-d3*inv_m+sqrt(s))/(1+inv_m^2)
          x5 <- (x3+y3*inv_m-d3*inv_m-sqrt(s))/(1+inv_m^2)
          y4 <- (d3+x3*inv_m+y3*inv_m^2+sqrt(s)*inv_m)/(1+inv_m^2)
          y5 <- (d3+x3*inv_m+y3*inv_m^2-sqrt(s)*inv_m)/(1+inv_m^2)

          junc$x[junc$atan == min(junc$atan)] <- x4
          junc$y[junc$atan == min(junc$atan)] <- y4
          junc$x[junc$atan == max(junc$atan)] <- x5
          junc$y[junc$atan == max(junc$atan)] <- y5

          junky <- rbind(junc, tmp_cell[tmp_cell$out == "out",])%>%
            mutate(id_cell = id_cell + 1,
                   my = mean(y),
                   mx = mean(x),
                   atan = atan2(y-my, x - mx))%>%
            plyr::arrange(atan)



          tmp_cello <- concavety(rbind(tmp_cell, junky)%>%mutate(id_point = paste0(round(x,8),";",round(y,8))))
        }else{

          #central point of the root hair apex
          x3 <- (x1+y1*m-d*m-sqrt(sigy))/(1+m^2)
          y3 <- (d+x1*m+y1*m^2-sqrt(sigy)*m)/(1+m^2)
          # perpendicular line
          inv_m <- -1/m
          d3 <- y3-inv_m*x3
          # distance
          s <- (hair_r^2)*(1+inv_m^2)-(y3-inv_m*x3-d3)^2
          #coordinate of the point
          x4 <- (x3+y3*inv_m-d3*inv_m+sqrt(s))/(1+inv_m^2)
          x5 <- (x3+y3*inv_m-d3*inv_m-sqrt(s))/(1+inv_m^2)
          y4 <- (d3+x3*inv_m+y3*inv_m^2+sqrt(s)*inv_m)/(1+inv_m^2)
          y5 <- (d3+x3*inv_m+y3*inv_m^2-sqrt(s)*inv_m)/(1+inv_m^2)

          junc$x[junc$atan == min(junc$atan)] <- x4
          junc$y[junc$atan == min(junc$atan)] <- y4
          junc$x[junc$atan == max(junc$atan)] <- x5
          junc$y[junc$atan == max(junc$atan)] <- y5

          junky <- rbind(junc, tmp_cell[tmp_cell$out == "out",])%>%
            mutate(id_cell = id_cell + 1,
                   my = mean(y),
                   mx = mean(x),
                   atan = atan2(y-my, x - mx))%>%
            arrange(atan)

          tmp_cello <- concavety(rbind(tmp_cell, junky)%>%mutate(id_point = paste0(round(x,8),";",round(y,8))))

        }

        tmp_cello$id_cell = tmp_cell$id_cell[1]
        # make a binding table with all root hair
        hairy <- rbind(hairy, tmp_cello)

      }else{
        id_h = id_h[id_h != h]
      }

    }


    rs1 <- rs1[rs1$id_cell %!in% id_h,]
    rs1$euc <- sqrt((x0-rs1$x)^2 + (y0-rs1$y)^2)
    rs1$out <- "none"
    rs8 <- rbind(rs1,hairy%>%select(-rank))

    return(rs8)
  }else{return(rs1)}


}
