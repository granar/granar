
#' @title Write .geo file for GMSH
#'
#' Write .geo file from cross section
#' geo file can be use in GMSH (Require GMSH https://gmsh.info/)
#' @param data The nodes dataframe of the cross section
#' @param path_geo The output file with .geo extension
#' @param Celldomain agrument to have a cell wall domain with filled cells (TRUE), or with hollow cells (FALSE)
#' @param dim The dimension of the input cross section, by default = 2
#' @keywords GMSH Geo
#' @export
#'
#'


write_geo <- function(data, dim = 2, path_geo, Celldomain =F){

  date = Sys.time()
  x1 = paste0('// Gmsh project created on ', date,'\nSetFactory("OpenCASCADE");\n//+\n')

  if(dim == 2){
    data$z = 0
    data$z1 = 0
    data$z2 = 0
  }

  k = h = j = 1
  txt = x1
  for(i in sort(unique(data$id_cell))){
    tmp = data%>%filter(id_cell == i)
    i_point = unique_point_id(tmp)

    i_point$idx = seq(k,-1+k+nrow(i_point),1)
    i_point$idx2 = c(i_point$idx[-1],i_point$idx[1])

    if(i < max(data$id_cell)){
      dots = paste0("Point(",i_point$idx,") = {",round(i_point$x,4), ' ,',round(i_point$y,4),' ,',round(i_point$z,4), ',1.0};\n', collapse = "")
      Lines = paste0("Line(",i_point$idx,") = {",i_point$idx, ", ", i_point$idx2,'};\n//+\n', collapse = "")
      Surf = paste0("Curve Loop(",j,") = {",paste0(i_point$idx, collapse = ", "),'};\n//+\nSurface(',
                    h,") = {", j,'};\n//+\n', collapse = "")
      if(Celldomain){
        Physical = paste0("Physical Surface(",h,") = {",h,"};\n")
      }else{
        Physical = paste0("//Physical Surface(",h,") = {",h,"};\n")
      }


      txt = paste0(txt, dots, Lines, Surf, Physical)
    }else{
      Physical_curve = paste0('Physical Curve("inner", 1) = {',paste0(1:(k-1), collapse = ", "),'};\n')
      dots = paste0("Point(",i_point$idx,") = {",round(i_point$x,4), ' ,',round(i_point$y,4),' ,',round(i_point$z,4), ',0.5};\n', collapse = "")
      Lines = paste0("Line(",i_point$idx,") = {",i_point$idx, ", ", i_point$idx2,'};\n//+\n', collapse = "")
      Surf = paste0("Curve Loop(",j,") = {",paste0(i_point$idx, collapse = ", "),'};\n//+\nPlane Surface(',
                    h,") = {", paste0(sort(seq(1,j,2), decreasing = T), collapse = ", "),'};\n//+\n', collapse = "")
      Physical = paste0("Physical Surface(0) = {",h,"};\n")

      txt = paste0(txt, dots, Lines, Surf, Physical, Physical_curve)
    }

    h = h+1
    j = j+2
    k = max(i_point$idx)+1
  }
  write(txt, path_geo)
}


#' @title remove duplicated points
#'
#'
#' @param tmp data to filter
#' @keywords GMSH Geo
#' @export
#'
#'

unique_point_id <- function(tmp){
  tmp$id_point = paste0('x:',tmp$x1,'_y:',tmp$y1,'_z:',tmp$z1)
  tmp$id_point2 = paste0('x:',tmp$x2,'_y:',tmp$y2,'_z:',tmp$z2)
  i_point = tibble(id_point = unique(c(tmp$id_point, tmp$id_point2)))

  x_coor = unlist(str_split(unlist(str_split(i_point$id_point, pattern = 'x:')),pattern = "_y"))
  y_coor = unlist(str_split(unlist(str_split(i_point$id_point, pattern = 'y:')),pattern = "_z"))
  z_coor = unlist(str_split(i_point$id_point, pattern = 'z:'))
  i_point$x = parse_number(x_coor[seq(2,length(x_coor),3)])
  i_point$y = parse_number(y_coor[seq(2,length(y_coor),3)])
  i_point$z = parse_number(z_coor[seq(2,length(z_coor),2)])

  i_point
  return(i_point)
}
