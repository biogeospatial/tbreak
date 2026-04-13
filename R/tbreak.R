library("Rbeast")
library("bfast")
library("zoo")
library("stringi")
library("lubridate")
library("methods")
library("parzer")
library("terra")
library("ncdf4")
library("sf")

get_date_vec_from_raster_names = function (raster) {
  #  generate time axis if needed
  #  assumes form 2024-01-24 somewhere in band name
  dates = time(raster)
  if (any (is.na(dates))) {
    dates = names(raster)
    pattern = r"(\b\d{4}-\d{2}-\d{2}\b)"
    m = regexpr(pattern, dates)
    dates = regmatches(dates, m)
    dates = strptime(dates, "%Y-%m-%d")
  }
  return(as.Date(dates))
}




assign_time_to_raster = function (raster, format= "%Y-%m-%d") {
  #  assumes form 2024-01-24 somewhere in band name
  #  dup from above - should be a function
  if (any (is.na(terra::time(raster)))) {
    dates = names(raster)
    if (any(stri_count_fixed(dates, "StdTime=") == 1)) {  #  ArcGIS mangled the dates
      base = lubridate::ymd("19000101", tz="UTC")
      pattern = r"(\d+\.\d+$)"
      m = regexpr(pattern, dates)
      dates = regmatches(dates, m)
      #  position of the 1 determines if it is a day or second - should use a day object
      dates = base - lubridate::ddays() + duration(as.numeric(dates), units="day")
    } else {  #  MODIS naming scheme
      pattern = r"(\b\d{4}-\d{2}-\d{2}\b)"
      m = regexpr(pattern, dates)
      dates = regmatches(dates, m)
      dates = strptime(dates, format, tz="UTC")
    }
    terra::time(raster) = dates
  }
  else {
    message ("Raster already has times")
  }
  return(raster)
}

get_cellnum_modis_coord = function (raster, coord) {
  coord = parse_coord_string(coord, crs(raster))

  cell_num = terra::cellFromXY (raster, cbind (x = coord[1], y = coord[2]))
  if (is.na(cell_num)) {
    stop ("Coord does not intersect the raster")
  }
  cell_num
}

get_ts_modis_coord = function (raster, coord) {
  cell_num = get_cellnum_modis_coord (raster, coord)
  #message (cell_num)
  unlist(raster[cell_num])
}

plot_ts_modis_coord = function (raster, coord, main=NULL) {
    # coord = parse_coord_string(coord)

    #  generate time axis if needed
    #  assumes form 2024-01-24 somewhere in band name
    #  dup from above - should be a function
    if (all (is.na(terra::time(raster)))) {
      dates = names(raster)
      pattern = r"(\b\d{4}-\d{2}-\d{2}\b)"
      m = regexpr(pattern, dates)
      dates = regmatches(dates, m)
      dates = strptime(dates, "%Y-%m-%d")
    }
    else {
      dates = terra::time(raster)
      #  nasty but we otherwise get bfastts errors
      dates = strptime(strftime(dates, format="%Y%m%d"), "%Y%m%d")
    }

    # cell_num = terra::cellFromXY (raster, cbind (x = coord[1], y = coord[2]))
    # if (is.na(cell_num)) {
    #   stop ("Coord does not intersect the raster")
    # }
    cell_num = get_cellnum_modis_coord (raster, coord)
    u = unlist(raster[cell_num])
    #u[u < -0.25] = NA  #  -0.3 is nodata, but not always

    if (all(is.na(u))) {
      message ("NA values only, skipping plot")
      return()
    }

    z = zoo::zoo(u, dates)
    plot (z, ylab = "Index", xlab = "Date", main = main)

    #  do we have any NAs?  highlight vals before and after
    na_prev_z = is.na(c(FALSE, z[-length(z)]))
    points (z[na_prev_z], col="red", pch = 16)
    na_next_z = c(is.na(z[-1]), FALSE)
    points (z[na_next_z], col="blue", pch=16)

    if (any(is.na(u))) {
      warning("NAs found in time series - break analyses might not work")
    }
    invisible(z)
}


point_to_cell_polygon = function (coord, raster) {
  coord = parse_coord_string(coord, crs(raster))

  c = res(raster)
  e = ext(raster)

  x1 = floor ((coord[1] - e$xmin) / c[1]) * c[1] + e$xmin
  x2 = x1 + c[1]
  y1 = floor ((coord[2] - e$ymin) / c[2]) * c[2] + e$ymin
  y2 = y1 + c[2]

  pol = sf::st_polygon(
    list(
      cbind(
        c(x1,x1,x2,x2,x1),
        c(y1,y2,y2,y1,y1)
      )
    )
  )
  pol = sf::st_sfc(pol, crs=crs(raster))
  pol
}

modis_crs = function () {
  m = r"(PROJCS["modis_sinusoidal",GEOGCS["GCS_Unknown_datum_based_upon_the_custom_spheroid",DATUM["D_Not_specified_based_on_custom_spheroid",SPHEROID["Custom_spheroid",6371007.181,0.0]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",0.0],UNIT["Meter",1.0]])"
  sf::st_crs(m)
}
