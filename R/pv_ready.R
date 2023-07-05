#' @title Format data for paraview visualisation
#'
#'
#' @param data node dataframe
#' @keywords root
#' @export
#' @examples
#' # rs <- pv_ready(rs1)
#'
#'

pv_ready <- function(data){

  data$neib1 = data$neib2 = data$neib3 <- NA
  data$neib1_type = data$neib2_type = data$neib3_type <- NA

  for(i in 1:nrow(data)){
    ne <- unique(data$id_cell[data$id_point == data$id_point[i] & data$id_cell != data$id_cell[i]])
    data$neib1[i] <- ne[1]
    data$neib1_type[i] <- as.character(data$type[data$id_cell == ne[1]][1])
    data$neib2[i] <- ne[2]
    data$neib2_type[i] <- as.character(data$type[data$id_cell == ne[2]][1])
    data$neib3[i] <- ne[3]
    data$neib3_type[i] <- as.character(data$type[data$id_cell == ne[3]][1])
  }

  junc_point <- data[!is.na(data$neib2) | (data$type == "epidermis" & !is.na(data$neib1)),]

  nodes <- data

  nodes<- nodes%>%
    dplyr::group_by(id_cell)%>%
    dplyr::filter(!duplicated(id_point))%>%
    dplyr::mutate(x1 = x,
                  y1 = y,
                  x2 = c(x[-1],x[1]),
                  y2 = c(y[-1],y[1]),
                  id_point2 = paste0(x2,";",y2),
                  w_length = sqrt((x1-x2)^2+(y1-y2)^2))

  on_go <- T
  while(on_go){

    last_point <- nodes%>%
      select(id_cell, last(starts_with("id_point")),
             last(starts_with("x")),
             last(starts_with("y"))) # take last colmn
    n <- colnames( last_point)
    n <- unique(readr::parse_number(n[-1])) # get the incremental value
    colnames(last_point) <- c("id_cell","id_point", "x","y") # generic col names

    last_point <- last_point%>%
      dplyr::group_by(id_cell)%>%
      dplyr::mutate(xx = ifelse(id_point %in% junc_point$id_point | id_point == "NA;NA",NA,c(x[-1],x[1])),
                    yy = ifelse(id_point %in% junc_point$id_point | id_point == "NA;NA",NA,c(y[-1],y[1])),
                    id_pointxy = paste0(xx,";",yy))%>%
      ungroup()%>%
      select(xx,yy,id_pointxy)

    on_go <- length(which(last_point$id_pointxy != "NA;NA"))>0
    c_name <- paste0(t(c("x", "y", "id_point")),n+1)
    colnames(last_point) <- c(c_name)

    nodes <- cbind(nodes%>%
                     ungroup(),last_point)
  }

  nodes <- nodes%>%
    dplyr::filter(id_point %in% junc_point$id_point,
                  w_length > 0)# only one wall per junction point

  # at this points all walls are double

  nodus <- NULL
  more <- nodes
  k <- 2
  while(nrow(more) > 0){


    last_x <- paste0("x",k) # the last point has coord_x
    last_y <- paste0("y",k) # the last point has coord_y
    tag_x <- paste0("x",k+1)  # the next possible points

    tmp <- more[is.na(more[,tag_x]),] # if next point is "nan" then
    tmp <- tmp%>%
      # dplyr::group_by(id_cell) %>%
      dplyr::mutate(xx = tmp[,last_x])%>%
      dplyr::mutate(yy = tmp[,last_y])

    deto <- tmp%>%select(-x1,-y1,-c(last_x),-c(last_y)) # remove points
    h <- 1
    # dealing with first and last point
    tmp_h <- tmp%>%
      mutate(xh = ifelse(x > xx, x, xx)) %>%
      mutate(yh = ifelse(x > xx, y, yy)) %>%
      mutate(yh = ifelse(x == xx, ifelse (y > yy, yy, y), yh))%>%
      mutate(xlh = ifelse(x > xx, xx, x)) %>%
      mutate(ylh = ifelse(x > xx, yy, y)) %>%
      mutate(ylh = ifelse(x == xx, ifelse (y > yy, y, yy), ylh))%>%
      select(xh,yh, xlh, ylh)
    colnames(tmp_h) <- c(paste0("x",h),paste0("y",h), paste0("x",k),paste0("y",k))
    deto <- cbind(deto, tmp_h)
    if(k >= 4){
      for(h in c(2:floor(k/2))){
        h_x <- paste0("x",h)
        lh_x <- paste0("x",k+1-h)
        h_y <- paste0("y",h)
        lh_y <- paste0("y",k+1-h)

        tmp_h <- tmp %>%
          mutate(hx = tmp[,h_x],
                 hy = tmp[,h_y],
                 lhx = tmp[,lh_x],
                 lhy = tmp[,lh_y])%>%
          mutate(xh = ifelse(x > xx, hx, lhx)) %>%
          mutate(yh = ifelse(x > xx, hy, lhy)) %>%
          mutate(xlh = ifelse(x > xx, lhx, hx)) %>%
          mutate(ylh = ifelse(x > xx, lhy, hy))%>%
          select(xh,yh,xlh,ylh)
        deto <- deto %>% select(-c(h_x,h_y,lh_x,lh_y))
        colnames(tmp_h) <- c(h_x,h_y,lh_x,lh_y)
        deto <- cbind(deto,tmp_h)
      }
    }
    nodus <- rbind(nodus,deto)
    k <- k + 1
    more <- more[!is.na(more[,tag_x]),] # reduce what left for the next loop
  }

  return(nodus)
}
