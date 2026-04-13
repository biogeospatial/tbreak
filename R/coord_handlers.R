library("stringi")

parse_coord_string = function (coord, out_crs=NULL) {
  coord = stri_trim(coord)
  
  in_metres = stri_sub(coord, -1) == "m"
  
  xy = numeric()
  #  "-33, 150"                       is nearmap lat/lon
  #  "-33,150"                        is QGIS lat/lon
  #  "150.7572609°E 33.9385950°S "    is ArcGIS lat/lon
  #  "13,907,112.15 -3,773,782.94 m " is ArcGIS projected
  #  "13907112.15,-3773782.94"        is QGIS projected
  #  '33°56'19.03"S 150°45'18.22"E'   is Google Earth lat/lon
  if (stri_count_regex(coord, ",\\s")) {
    xy = stri_split_regex(coord, ",\\s")[[1]][1:2]
    x = as.numeric (xy[1])
    if (x > 360 || x < -180) {
      in_metres = TRUE
    } else {
      xy[1:2] = xy[2:1]  #  lat/long to x/y
    }
  } else if (stri_count_regex(coord, "\\s")) {
    xy = stri_split_regex(coord, "\\s")[[1]][1:2]
    xy = stri_replace_all_fixed (xy, replacement="", pattern=",")
    if (stri_count_regex(xy[2], "[Ee]")) {
      xy[1:2] = xy[2:1]  #  convert to lon/lat
    }
  } else {
    xy = stri_split_fixed(coord, ",")[[1]][1:2]
    x = as.numeric (xy[1])
    if (x > 360 || x < -180) {
      in_metres = TRUE
    } else {
      xy[1:2] = xy[2:1]  #  lat/long to x/y
    }
  }
  
  if (!in_metres) {
    library("parzer")
    library("sf")
    
    poss_plain_dd = stri_detect (coord, regex = r"(^-?\d+\.\d+,\s*-?\d+\.\d+$)")
    
    if (poss_plain_dd) {
      xy = stri_split_regex(coord, ",\\s+")[[1]][2:1]
    }
    x = parse_lon(xy[1])
    y = parse_lat(xy[2])
    point = st_sfc(st_point(c(x,y)), crs = 4326)
    if (!(is.null(out_crs) || out_crs=="")) {
      point = st_transform (point, out_crs)
    }
    p = st_coordinates(point)
    return (c(p[1], p[2]))
  }
  
  x_hemi = stri_extract(xy[1], regex = "[NESW]")
  if (is.na(x_hemi)) {
    x_hemi = "E"
  }
  x = stri_replace_all_fixed(xy[1], replacement = "", pattern = x_hemi)
  x = as.numeric(x)
  if (x_hemi == "W") {
    x = x * -1
  }
  
  y_hemi = stri_extract(xy[2], regex = "[NESW]")
  if (is.na(y_hemi)) {
    y_hemi = "N"
  }
  y = stri_replace_all_fixed(xy[2], replacement = "", pattern = y_hemi)
  y = as.numeric(y)
  if (y_hemi == "S") {
    y = y * -1
  }
  
  return (c(x, y))
}


#  convert a coordinate copied from ArcGIS to a format usable with terra
#  requires the extent be stored on the rbeast object
coord2idx_rbeast = function (b, coord) {
  coord = parse_coord_string(coord, b$crs)
  x = coord[1]
  y = coord[2]
  
  extent = terra::ext(b$ext)
  
  if (x < extent$xmin || x >= extent$xmax) {
    stop ("X coord is outside data set extent")
  }
  if (y < extent$ymin || y >= extent$ymax) {
    stop ("Y coord is outside data set extent")
  }
  
  ncol = b$ncols
  nrow = b$nrows
  
  c = 1    + floor ((x - extent$xmin) / (extent$xmax - extent$xmin) * ncol)
  r = nrow - floor ((y - extent$ymin) / (extent$ymax - extent$ymin) * nrow)
  
  res = c(r, c, use.names=FALSE)
  return (res)
}

