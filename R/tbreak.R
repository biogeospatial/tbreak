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

get_date_vec_from_raster = function (raster) {
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
  return(dates)
}

calc_and_plot_beast_modis_coord = function (raster, coord, main = NULL, start_time=NULL, ...) {
  coord = parse_coord_string(coord)

  dates = get_date_vec_from_raster(raster)

  cell_num = terra::cellFromXY (raster, cbind (x = coord[1], y = coord[2]))
  if (is.na(cell_num)) {
    stop ("Coord does not intersect the raster")
  }
  Y = unlist(raster[cell_num])
  Y[Y < -0.25] = NA

  metadata = list(
    time             = dates,
    isRegularOrdered = FALSE,    # IRREGULAR input
    #whichDimIsTime   = 3,        # 437 is the ts length, so set it to '3' here.
    # time$datestr     = datestr,  # date info is contained in the file names
    # time$strfmt      = 'LT05_018032_20080311.yyyy-mm-dd',
    deltaTime        = 16/365,     # MODIS data are 16 days
    #deltaTime        = 1/12,
    #    period = 32/365
    period = 1
    #period           = 16/365
    #startTime        = ifelse (is.null(start_time, time(tmp_ras)[1], as.Date(start_time)))
  )

  #  minimal for now
  extra = list (
    numThreadsPerCPU = 3,
    numParThreads    = 30
  )

  o = Rbeast::beast123 (Y, metadata = metadata, extra = extra, ...)

  if (!is.na(o$R2)) {
    plot(o, main = main)
  } else {
    message ("Unable to fit model, no plot available")
  }

  return (o)
}

tiled_beast_modis = function (raster, tile_size=64, start_time = NULL, ...) {
  res = res(raster)[1:2] * tile_size
  ext = ext(raster)
  nrows = ceiling((ext[2] - ext[1]) / res[2])
  ncols = ceiling((ext[4] - ext[3]) / res[1])
  #  way too many args but terra does odd things at the moment
  #  and this does what we want
  rr = rast(
    nrows = nrows,
    ncols = ncols,
    xmin  = ext[1],
    xmax  = ext[2],
    ymin  = ext[3],
    ymax  = ext[4],
    crs   = crs(raster),
    res   = res(raster)[1:2] * tile_size
  )
  v = as.polygons(rr)
  v = terra::intersect(v, ext(raster))
  v$beast_id = 1:nrow(v)
  v = st_as_sf(v)  #  make it an sf object

  rm (rr)
  gc()

  b = list()
  subset = 1:nrow(v)
  #subset = 1:3  #  for debug
  for (i in v$beast_id[subset]) {
    if (is.na(i)) {break}  #  for debug
    r = crop(raster, v[i,], ext=TRUE)
    b[[i]] = beast_modis(r, start_time=start_time, ...)
  }

  #  maybe convert v to an sf object?
  bm = list (index = v, beasts = b)

  return (bm)
}

beast_modis = function (raster, start_time = NULL, ...) {

  #  work on a copy in case it is not in memory
  #  then the copy should be GC'd when that is called
  tmp_ras = raster
  tmp_ras = toMemory(tmp_ras)   # To use beast, make sure all the data is read into memory
  dims    = dim(tmp_ras)

  dates = get_date_vec_from_raster(raster)

  # Y = values(tmp_ras)
  #dim(Y)   = dims[c(2,1,3)]    # Assign column-major dim expected by Rbeast
  #  avoid the need to transpose
  Y = as.matrix(tmp_ras, wide=TRUE)
  dim(Y)   = dims[c(1,2,3)]  #  could just use dims directly...

  metadata = list(
    time             = dates,
    isRegularOrdered = FALSE,    # IRREGULAR input
    whichDimIsTime   = 3,        # 437 is the ts length, so set it to '3' here.
    # time$datestr     = datestr,  # date info is contained in the file names
    # time$strfmt      = 'LT05_018032_20080311.yyyy-mm-dd',
    deltaTime        = 16/365,     # MODIS data are 16 days
    #deltaTime        = 1/12,
    #    period = 32/365
    period = 1
    #period           = 16/365
    #startTime        = ifelse (is.null(start_time, time(tmp_ras)[1], as.Date(start_time)))
  )

  #  minimal for now
  extra = list (
    numThreadsPerCPU = 3,
    numParThreads    = 30
  )

  #browser()
  tmp_ras = NULL
  gc()


  o = Rbeast::beast123 (Y, metadata = metadata, extra = extra, ...)

  #  pack in some extra info - underhanded but works for now
  o$ext = as.vector(ext(raster)) #  ext objects do not survive serialisation otherwise
  o$crs = crs(raster)

  gc()
  return (o)
}

parse_coord_string = function (coord) {
  coord = stri_trim(coord)

  in_metres = stri_sub(coord, -1) == "m"

  xy = stri_split_regex(coord, "\\s")[[1]][1:2]
  xy = stri_replace_all_fixed (xy, replacement="", pattern=",")

  if (!in_metres) {
    library("parzer")
    library("sf")

    poss_plain_dd = stri_detect (coord, regex = r"(^-?\d+\.\d+,\s+-?\d+\.\d+$)")

    if (poss_plain_dd) {
      xy = stri_split_regex(coord, ",\\s+")[[1]][2:1]
    } else {
      x_hemi = toupper (stri_extract(xy[1], regex = "[NESWnesw]"))
      if (x_hemi %in% c("N", "S")) {
        xy = rev (xy)
      }
    }
    x = parse_lon(xy[1])
    y = parse_lat(xy[2])
    point = st_sfc(st_point(c(x,y)), crs = 4326)
    point2 = st_transform (point, modis_crs())
    p = st_coordinates(point2)
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


ext_from_arcgis_coord = function (coord, xoff=100000, yoff=-100000) {
  coord = parse_coord_string(coord)
  x1 = coord[1]
  x2 = coord[1] + xoff
  y1 = coord[2]
  y2 = coord[2] + yoff
  terra::ext(min(x1, x2), max(x1, x2), min(y1, y2), max(y1, y2))
}

#  convert a coordinate copied from ArcGIS to a format usable with terra
#  requires the extent be stored on the rbeast object
coord2idx_rbeast = function (b, coord) {
  coord = parse_coord_string(coord)
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


plot_beast_modis_coord = function (b, coord, t=FALSE, main=NULL) {

  if (isa_tiled_beast(b)) {
    p = st_point (parse_coord_string(coord))
    target_tile = st_intersects(p, b$index)[[1]]
    return (plot_beast_modis_coord(b$beasts[[target_tile]], coord, t, main))
  }

  rowcol = coord2idx_rbeast(b, coord)


  #  transpose for beast coords
  if (t) {
    rc = rowcol[c(2,1)]
  }
  else {
    rc = rowcol[c(1,2)]
  }
  message ("Plotting image row ", rc[1], ", col ", rc[2])
  message ("Mean number of trend change points is ",  b$trend$ncp[rc[1],rc[2]])
  message ("Mean number of season change points is ", b$season$ncp[rc[1],rc[2]])
  plot (b[rc], main=main)

}

load_data = function (file = NULL, drivers=NULL) {
  if (is.null(file)) {
    file = file.choose()
  }
  r = rast(file, drivers=drivers)
  if (any (is.na(terra::time(r)))) {
    dates = strptime(names(r), "%Y-%m-%d")
    if (any(is.na(dates))) {
      stop ("Some of the field names do not satisfy the date format requirement (yyyy-mm-dd)")
    }
    terra::time(r) = as.Date(dates)
  }
  return (r)
}


beast_time_to_date = function (x) {
  #  this fails because the current month and day are assumed for %Y on its own
  #b$time_as_date = as.Date(as.character(floor(b$time)), format="%Y") + (b$time - floor (b$time))*365
  as.Date(paste(floor(x), ceiling((x - floor (x))*365), sep=""), format="%Y%j")
}

#  need to do one part at a time
rasterise_tiled_beast = function (tb) {
  if (is.null(tb$index)) {
    stop ("An index must be set on the tiled beast object")
  }

  #b$time_as_date = lubridate::date_decimal (b$time)

  results = list ()

  max_idx = max(tb$index$beast_id)

  for (idx in tb$index$beast_id) {
    b = tb$beasts[[idx]]
    for (component in c("trend", "season")) {
      for (subcomponent in names(b[[component]])) {
        if (!is.null (b[[component]][[subcomponent]])) {
          message (sprintf("tile %s of %s: %s: %s", idx, max_idx, component, subcomponent))
          label = sprintf ("%s_%s", component, subcomponent)
          results[[label]][[idx]] = beastbit2raster(b, component, subcomponent)
        }
      }
    }
    for (component in c("R2", "RMSE", "sig2", "marg_lik")) {
      message (sprintf("tile %s of %s: %s", idx, max_idx, component))
      label = component
      results[[label]][[idx]] = beastbit2raster(b, component)
    }
  }

  #  now make mosaics from spat raster collections
  for (name in names(results)) {
    results[[name]] = terra::mosaic(terra::sprc(results[[name]]))
  }

  results
}

rasterise_beast = function (b) {
  if (is.null(b$ext)) {
    stop ("An extent must be set on the beast object")
  }

  b$time_as_date = lubridate::date_decimal (b$time)

  results = list ()
  #trend_ncp = beastbit2raster (b, "trend", "ncp")

  for (component in c("trend", "season")) {
    for (subcomponent in names(b[[component]])) {
      if (!is.null (b[[component]][[subcomponent]])) {
        message (sprintf("%s: %s", component, subcomponent))
        label = sprintf ("%s_%s", component, subcomponent)
        results[[label]] = beastbit2raster(b, component, subcomponent)
      }
    }
  }
  for (component in c("R2", "RMSE", "sig2", "marg_lik")) {
    message (component)
    label = component
    results[[label]] = beastbit2raster(b, component)
  }


  results
}

get_beast_component_list = function (b) {

  results = list()

  if (isa_tiled_beast(b)) {
    #  tiled beast
    #  loop over the lot because one day
    #  we might store nulls for nodata blocks,
    #  but break once we have found something
    for (idx in b$index$beast_id) {
      bb = b$beasts[[idx]]
      r = get_beast_component_list(bb)
      results = c(results, r)
      if (length(names(results)) > 0) {
        break
      }
    }
    return (results)
  }

  for (component in c("trend", "season")) {
    for (subcomponent in names(b[[component]])) {
      if (!is.null (b[[component]][[subcomponent]])) {
        results[[component]][[subcomponent]] = TRUE
      }
    }
  }
  for (component in c("R2", "RMSE", "sig2", "marg_lik")) {
    results[[component]] = TRUE
  }


  results
}

beastbit2raster = function (b, component = "trend", subcomponent = "ncp", template=NULL) {

  if (isa_tiled_beast(b)) {
    rasters = list()

    for (idx in b$index$beast_id) {
      rr = beastbit2raster(b$beasts[[idx]], component, subcomponent, template)
      if (max(minmax(rr)) == Inf) {  #  terra::max does not exist at the moment
        rr[rr == Inf] = NA
      }
      rasters[[idx]] = rr
    }
    return(mosaic(sprc(rasters)))
  }

  if (class(b) != "beast") {
    stop ("Need an Rbeast result")
  }

  temporal_subcomponents = c(
    "cpOccPr", "Y", "SD", "order",
    "slp",     "slpSD",
    "slpSgnPosPr", "slpSgnZeroPr",
    "amp",     "ampSD"
  )

  valid_components = c("trend", "season", "R2", "RMSE", "sig2", "marg_lik")
  if (!(component %in% valid_components)) {
    stop (cat ("component args must be in list: ", paste(valid_components, sep=", ")))
  }
  if ((component %in% c("trend", "season")) && is.null(b[[component]][[subcomponent]])) {
    stop (sprintf ("Cannot find b$%s$%s, or it is null", component, subcomponent))
  }

  ext = b$ext
  if (is.null(ext) && (is.null (template) || !hasMethod("ext", class(template)))) {
    stop ("No extent found. Need an Rbeast object with an ext slot, or a template object with an ext method")
  }
  ext = ext(ext) #  ensure we have an extent object

  crs = b$crs
  if (is.null(crs) && (is.null (template) || !hasMethod("crs", class(template)))) {
    stop ("No coord sys found.  Need an Rbeast object with a crs slot, or a template object with a crs method")
  }

  rasterise = function (data) {
    m = matrix(
      data,
      b$nrows[1],
      b$ncols[1],
      byrow = FALSE
    )
    r = terra::rast(m, extent=ext, crs=crs)
  }

  r = NULL

  if (component %in% c("trend", "season")) {
    subdata = b[[component]][[subcomponent]]
    dims = dim(subdata)
    if (length (dims) > 4) {
      message (sprintf ("Too many dims, skipping: %s, %s", component, subcomponent))
    }
    else if (length(dims) < 3) {
      r = rasterise (subdata)
      names(r) = sprintf ("%s_%s", component, subcomponent)
    }
    else if (length(dims) == 3) {
      nbands = dims[3]
      r = do.call("c", lapply(1:nbands, function(i){
        return (rasterise (subdata[,,i]))
      }))

      #  These should be temporal data.
      if (subcomponent %in% temporal_subcomponents && nbands == length(b$time)) {
        #message (sprintf ("Temporal data for %s, %s", component, subcomponent))
        t = b$time_as_date  #  use pre-calculated if exists
        if (is.null(t)) {
          t = lubridate::date_decimal (b$time)
        }
        terra::time(r) = t
      }

      names(r) = paste0 (
        sprintf ("%s_%s", component, subcomponent),
        formatC(1:nbands, flag="0", width=nchar(nbands))
      )

    }
    else {
      ci_names = c("lower", "upper")
      r = list()
      for (j in 1:2) {

        nbands = dims[3]
        rr = do.call("c", lapply(1:nbands, function(i){
          return (rasterise (subdata[,,i,j]))
        }))

        names(rr) = paste0 (
          sprintf ("%s_%s_%s", component, subcomponent, ci_names[j]),
          formatC(1:nbands, flag="0", width=nchar(nbands))
        )
        if (length(r) > 0) {
          r = append (r, rr)
        }
        else {
          r = rr
        }
      }

    }
  }
  else {
    subdata = b[[component]]
    dims = dim(subdata)
    r = rasterise (subdata)
    names(r) = component
  }

  r
}

#  pretty basic check until we develop a tiled beast class
isa_tiled_beast = function(b) {
  if (is.atomic(b) || (is.null(b$index) && is.null(b$beasts))) {
    return (FALSE)
  }
  return (all(sapply (b$beasts, FUN=function(x){is.null(x) || class(x)=='beast'})))
}


export_beast_rasters = function (b, dir, prefix="", overwrite=FALSE) {

  name_list = get_beast_component_list(b)
  message ("Exporting now")

  outputs = character()

  for (component_name in names(name_list)) {
    subcomponent = name_list[[component_name]]
    if (is.list(subcomponent)) {
      for (subcomponent_name in names(subcomponent)) {
        name = sprintf ("%s_%s", component_name, subcomponent_name)
        pfx = file.path (dir, paste0(prefix, name))
        outfile = paste0 (pfx, ".tif")
        message (outfile)

        r = beastbit2raster(b, component_name, subcomponent_name)
        if (max(minmax(raster)) == Inf) {
          r[r == Inf] = NA
        }

        writeRaster(r, outfile, overwrite=overwrite)
        outputs = c(outputs, outfile)

        #  also dump netCDF for temporal data
        if (!anyNA(terra::time(r))) {
          outfile = paste0 (pfx, ".nc")
          message (outfile)
          writeCDF(r, outfile, overwrite=overwrite)
          outputs = c(outputs, outfile)
        }
      }
    }
    else {
      name = sprintf ("%s", component_name)
      pfx = file.path (dir, paste0(prefix, name))
      outfile = paste0 (pfx, ".tif")
      message (outfile)

      r = beastbit2raster(b, component_name)
      if (max(minmax(raster)) == Inf) {
        r[r == Inf] = NA
      }

      writeRaster(r, outfile, overwrite=overwrite)
      outputs = c(outputs, outfile)
    }
  }

  invisible(outputs)
}



plot_bfast_modis_coord = function (raster, coord, h=0.15, main=NULL) {
  coord = parse_coord_string(coord)

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

  cell_num = terra::cellFromXY (raster, cbind (x = coord[1], y = coord[2]))
  if (is.na(cell_num)) {
    stop ("Coord does not intersect the raster")
  }
  u = unlist(raster[cell_num])
  u[u < -0.25] = NA

  if (sum(is.na(u)) > (length(u) / 10)) {
    message ("More than 10% of records are NA, skipping bfast generation")
    return()
  }

  t2 = bfast::bfastts(u, dates, type = '16-day')
  tb = bfast::bfast(t2, h=h)
  plot(tb, main=main)
  invisible (tb)
}

assign_time_to_raster = function (raster, format= "%Y-%m-%d") {
  #  assumes form 2024-01-24 somewhere in band name
  #  dup from above - should be a function
  if (any (is.na(terra::time(raster)))) {
    dates = names(raster)
    pattern = r"(\b\d{4}-\d{2}-\d{2}\b)"
    m = regexpr(pattern, dates)
    dates = regmatches(dates, m)
    dates = strptime(dates, format)
    terra::time(raster) = dates
  }
  else {
    message ("Raster already has times")
  }
  return(raster)
}


plot_ts_modis_coord = function (raster, coord, main=NULL) {
    coord = parse_coord_string(coord)

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

    cell_num = terra::cellFromXY (raster, cbind (x = coord[1], y = coord[2]))
    if (is.na(cell_num)) {
      stop ("Coord does not intersect the raster")
    }
    u = unlist(raster[cell_num])
    u[u < -0.25] = NA  #  -0.3 is nodata

    if (all(is.na(u))) {
      message ("NA values only, skipping plot")
      return()
    }

    z = zoo::zoo(u, dates)
    plot (z, xlab = "Index", ylab = "Date", main = main)

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
  coord = parse_coord_string(coord)

  c = res(raster)
  e = ext(raster)

  x1 = floor ((coord[1] - e$xmin) / c[1]) * c[1] + e$xmin
  x2 = x1 + c[1]
  y1 = floor ((coord[2] - e$ymin) / c[2]) * c[2] + e$ymin
  y2 = y1 + c[2]

  pol = st_polygon(
    list(
      cbind(
        c(x1,x1,x2,x2,x1),
        c(y1,y2,y2,y1,y1)
      )
    )
  )
  pol = st_sfc(pol, crs=modis_crs())
  pol
}

modis_crs = function () {
  m = r"(PROJCS["modis_sinusoidal",GEOGCS["GCS_Unknown_datum_based_upon_the_custom_spheroid",DATUM["D_Not_specified_based_on_custom_spheroid",SPHEROID["Custom_spheroid",6371007.181,0.0]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",0.0],UNIT["Meter",1.0]])"
  st_crs(m)
}
