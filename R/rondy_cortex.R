#' @title Add inter cellular space
#'
#' Re-do the cortex layer and make it more round shaped
#' @param params The input dataframe
#' @param all_cells The cellular dataframe
#' @param center The cross-section center
#' @keywords root
#' @export
#' @examples
#' # all_cells <- rondy_cortex(params, all_cells, center)
#'

rondy_cortex <- function(params, all_cells, center){

  if( is.null(params$value[params$name == "inter_cellular_space" & params$type == "size"])){
    warning("Please specify the inter cellular space component in the input file")
    return(NULL)
  }

  random_fact <- params$value[params$name == "randomness"] / 10
  cor_d <- params$value[params$name == "cortex" & params$type == "cell_diameter"]

  all_cortex <- all_cells[all_cells$type %in% c("cortex","endodermis", "exodermis"),]

  icp_size <- params$value[params$name == "inter_cellular_space" & params$type == "size"]
  if(length(icp_size)> 0){
    scaling <- 1-(icp_size/cor_d)
    if(scaling >= 0.99){
      scaling = 0.99
    }
  }else{scaling <- 0.95}


  if(length(all_cells$id_group[all_cells$type == "xylem"]) > 0){
    k_max_xylem <- max(all_cells$id_group[all_cells$type == "xylem"])
  }else{k_max_xylem <- 0}

  ctess <- deldir(all_cortex$x, all_cortex$y, digits = 8)
  idc <- unique(all_cortex$id_cell)
  idc <- 1:length(idc)
  rc <- ctess$dirsgs[ctess$dirsgs$ind1 %in% idc |
                       ctess$dirsgs$ind2 %in% idc,]
  rc <- rc%>% arrange(ind1)
  rc2 <- data.frame(x = rc$x1, y=rc$y1, id_cell = rc$ind1)
  rc2 <- rbind(rc2, data.frame(x = rc$x2, y=rc$y2, id_cell = rc$ind1))
  rc2 <- rbind(rc2, data.frame(x = rc$x2, y=rc$y2, id_cell = rc$ind2))
  rc2 <- rbind(rc2, data.frame(x = rc$x1, y=rc$y1, id_cell = rc$ind2))

  inner <- min(all_cortex$radius[all_cortex$type %in% c("cortex")]) # [all_cortex$type %in% c("endodermis", "cortex")]
  outer <- max(all_cortex$radius[all_cortex$type %in% c("cortex")]) # no intercellular space between exo cortex and cortex

  rc2 <- rc2%>%mutate(euc = sqrt((x-center)^2+(y-center)^2))%>%
    dplyr::group_by(id_cell)%>%
    dplyr::mutate(mx = mean(x),
                  my = mean(y),
                  atan = atan2(y-my, x - mx))%>%
    dplyr::arrange(id_cell, atan)

  all_cortex$id_cell <- 1:nrow(all_cortex)
  rc1 <- merge(rc2, all_cortex[,c("id_cell", "type", "radius", "id_layer")], by="id_cell")

  rcin <- rc2%>%
    filter(euc > inner,
           euc < outer)%>%
    mutate(ID = paste0(x,y))%>%
    filter(!duplicated(ID))

  all_inter <- data.frame(angle = ifelse(rcin$y-center >= 0, acos((rcin$x - center)/rcin$euc),
                                         2*pi-acos((rcin$x - center)/rcin$euc)),
                          radius = rcin$euc,
                          x = rcin$x, y = rcin$y,
                          id_layer = all_cortex$id_layer[1]+0.5,
                          id_cell = 1:nrow(rcin),
                          type = "inter_cellular_space",
                          order = params$value[params$name == "cortex" & params$type == "order"]+0.5,
                          id_group = 0 # if too close, they should be merge but not now
  )
  if(length(params$value[params$name =="inter_cellular_space"]) > 0){
    coef_icp <- (10 * params$value[params$name =="inter_cellular_space" & params$type == "ratio"]) # here to modulate icp proportion coeficient
    if(coef_icp > 1){coef_icp = 1}
    inter_cellular_proportion <- coef_icp*nrow(all_inter)
  }else{
    inter_cellular_proportion <- 0.5*nrow(all_inter)
  }
  if(inter_cellular_proportion == 0){
    all_inter <- NULL
  }else{
    to_keep <- sample(1:nrow(all_inter), round(inter_cellular_proportion), replace=F)
    all_inter <- all_inter[all_inter$id_cell %in% to_keep,]
    all_inter$id_point <- paste0(all_inter$x,";",all_inter$y)
    all_inter <- all_inter%>%
      filter(!duplicated(id_point))%>%
      select(-id_point)
  }

  nodes <- vertex(rc1%>%filter(type == "cortex"))
  nodes <- nodes %>%
    filter(wall_length > 0)%>%
    mutate(m = (y2-y1)/(x2-x1),
           k = y1-m*x1,
           r_dist = abs(k+m*mx-my)/sqrt(1+m^2)) # distance between a point (mx,my) to a segment defined by to point (x1,y1; x2,y2)

  cor <- nodes%>%
    filter(wall_length > 0)%>%
    dplyr::group_by(id_cell)%>%
    dplyr::mutate(radius = min(r_dist))%>%
    filter(!duplicated(id_cell))

  cor$radius[cor$id_layer %in% c(inner, outer)] <- cor$radius[cor$id_layer %in% c(inner, outer)]*0.45

  cor$id_group = 1:nrow(cor)

  circus <- seq(-0.95,0.95,0.95/4)
  cir <- data.frame(x_cir = rep(circus,2*nrow(cor)))%>%
    mutate(y_cir = rep(c(sqrt(1-circus^2),-sqrt(1-circus^2)),nrow(cor)),
           id_group = sort(rep(1:nrow(cor),2*length(circus))))

  cor_frontier <- merge(cir, cor[,c("id_group", "radius", "mx", "my", "id_layer")], by = "id_group")%>%
    transmute(radius = radius,
              x = x_cir*radius*scaling+mx,
              y = y_cir*radius*scaling+my,
              euc = sqrt((mx-center)^2+(my-center)^2),
              angle = ifelse(my-center > 0,acos((mx - center)/euc),
                             2*pi-acos((mx - center)/euc)) ,
              id_layer = id_layer,
              id_cell = 1,
              type = "cortex",
              order = params$value[params$name == "cortex" & params$type == "order"],
              id_group = id_group
    )%>%
    select(-euc)

  all_cells <- rbind(all_cells[all_cells$type != "cortex",], cor_frontier)
  all_cells$id_group[all_cells$type == "cortex"] <- all_cells$id_group[all_cells$type == "cortex" & all_cells$id_group != 0] + k_max_xylem
  all_cells <- rbind(all_cells, all_inter)

  # reset the cell ids
  all_cells$id_cell <- c(1:nrow(all_cells))
  return(all_cells)
}
