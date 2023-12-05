#' @title Do circle packing to add Xylem element in stele
#'
#'
#' @param params input parameter
#' @param center center of the root cross section
#' @param all_cells dataframe with the id_group of the cells
#' @keywords root
#' @export
#' @examples
#' # rs <- septa(rs1)
#'
#'

pack_xylem <- function(all_cells, params, center){


  if(length(all_cells$id_group[all_cells$type == "cortex"])> 0){
    k_max_cortex <- max(all_cells$id_group[all_cells$type == "cortex"])
  }else{k_max_cortex  <- 0}

  nX <- params$value[params$name == "xylem" & params$type == "n_cells"]
  xylarea <- nX* pi * (params$value[params$name == "xylem" & params$type == "cell_diameter"]/2)^2

  # stele diameter = before pericycle
  stele_diameter<-params$value[params$name == "stele" & params$type == "layer_diameter"]

  # Xylem cell size -> beta distribution
  stele_cell_area <- pi*(params$value[params$name == "stele" & params$type == "cell_diameter"]/2)^2
  max_Xylem_cell_area = pi*(params$value[params$name == "xylem" & params$type == "max_size"]/2)^2

  # Nparenchyma = stele_area -
  nparenchyma<-(((((stele_diameter/2))^2)*pi)-xylarea)/stele_cell_area

  # Correction due to bugs if too much parenchyma
  if(nparenchyma>=1000){
    nparenchyma=nparenchyma*0.86}
  if(nparenchyma>=300 & nparenchyma<1000){
    nparenchyma=nparenchyma*0.95}

  # Define stell cell diameter SD
  if(length(params$value[params$name == "stele" & params$type == "SD"])==0){
    stele_cell_area_SD <- 0.1*stele_cell_area*10^6
  }else{
    stele_cell_area_SD <- params$value[params$name == "stele" & params$type == "SD"]
  }

  # Areas of parenchyma cells and xylem
  areas <- c(rnorm(round(nparenchyma), stele_cell_area*10^6, stele_cell_area_SD),
             rbeta(nX,shape1 = 2, shape2 = 4)*
               max_Xylem_cell_area*10^6) #rnorm (number, mean, SD) rbeta (number, shape1, shape2)* area max
  areas<-as.data.frame(areas)
  areas$type="NA"
  areas$type[1:round(nparenchyma)]="parenchyma"
  areas$type[round(nparenchyma+1):nrow(areas)]="xylem"

  areas<-rbind(areas%>%filter(type == "xylem"),areas%>%filter(type != "xylem"))# areas_phlo

  # Bind xylem and parenchyma in a way that xylem is in the 10% first cells.
  # Since first cells will be placed in the center, that allows xylem to be placed in the 10% central cells
  areas <- rbind(areas[sample(1:round(0.1*nrow(areas))),],areas[c((round(0.1*nrow(areas))+1):nrow(areas)),])

  # PACK CIRCLES -> arrange tangentialy the (circle) cells with no overlapping
  # the first cells of the Areas are arranged first from the center.
  areas$areas <- abs(areas$areas)
  areas$id <- 1:nrow(areas)
  # Generate the layout
  packing <- as.data.frame(circleProgressiveLayout(areas,sizecol ="areas"))
  packing <- cbind(packing,areas)
  dat.gg <- circleLayoutVertices(packing, npoints=25, sizetype = "radius")
  dat.gg<-full_join(dat.gg,areas,by="id")


  # Store packed cells
  packing <- as.data.frame(packing)%>%
    mutate(radius = radius/1000,
           x = x/1000 + center,
           mx = x,
           y= y /1000+ center,
           my = y,
           area = pi*radius^2,
           type = type,
           id_cell = id)

  # ======================================================================================================
  # Make xylem frontiers
  xyl <- packing[packing$type == "xylem", ]

  xyl$id_group = 1:nrow(xyl)
  xyl$id_layer = 1.5
  circus <- seq(-0.95,0.95,0.95/4)
  cir <- data.frame(x_cir = rep(circus,2*nrow(xyl)))%>%
    mutate(y_cir = rep(c(sqrt(1-circus^2),-sqrt(1-circus^2)),nrow(xyl)),
           # mx = rep(cor$mx,2*length(circus)),
           # my = rep(cor$my,2*length(circus)),
           id_group = sort(rep(1:nrow(xyl),2*length(circus))))

  scaling = 0.95

  xyl_frontier <- merge(cir, xyl[,c("id_group", "radius", "mx", "my", "id_layer")], by = "id_group")%>%
    transmute(radius = radius,
              x = x_cir*radius*scaling+mx, # Check scaling in rondicortex (0.95)
              y = y_cir*radius*scaling+my,
              euc = sqrt((mx-center)^2+(my-center)^2),
              angle = ifelse(my-center > 0,acos((mx - center)/euc),
                             2*pi-acos((mx - center)/euc)) ,
              id_layer = id_layer,
              id_cell = 1,
              type = "xylem",
              order = params$value[params$name == "xylem" & params$type == "order"],
              id_group = id_group
    )%>%
    select(-euc)


  # ======================================================================================================
  # filter without xyl and check column names
  packing_par <- subset(packing, packing$type =="parenchyma" |packing$type =="phloem"|packing$type =="test") %>%
    transmute(radius = radius,
              x=x,
              y=y,
              euc = sqrt((x-center)^2+(y-center)^2),
              angle = ifelse(y-center > 0,acos((x - center)/euc),
                             2*pi-acos((x - center)/euc)),
              id_layer = 1.5,
              id_cell = 1,
              type = type,
              order = 1,
              id_group = 0
    )%>%
    filter(euc < stele_diameter/2-params$value[params$name == "stele" & params$type == "cell_diameter"]/2)%>%
    select(-euc)

  # ======================================================================================================
  packing_par$id_cell<-1:length(packing_par$id_cell)
  xyl_frontier$id_cell<-(length(packing_par$id_cell)+1):(length(packing_par$id_cell)+length(xyl_frontier$id_cell))
  packing <- rbind (packing_par , xyl_frontier)

  packing$id_group[packing$type == "xylem"] <- packing$id_group[packing$type == "xylem" & packing$id_group != 0] + k_max_cortex
  return(packing)
}
