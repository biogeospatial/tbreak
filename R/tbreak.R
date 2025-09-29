beast_modis = function (raster, start_time = NULL, ...) {
  
  #  work on a copy in case it is not in memory
  #  then the copy should be GC'd when that is called
  tmp_ras = raster
  tmp_ras = toMemory(tmp_ras)   # To use beast, make sure all the data is read into memory
  dims    = dim(tmp_ras)
  
  #  generate time axis if needed
  #  assumes form 2024-01-24 somewhere in band name
  if (all (is.na(time(tmp_ras)))) {
    dates = names(tmp_ras)
    pattern = r"(\b\d{4}-\d{2}-\d{2}\b)"
    m = regexpr(pattern, dates)
    dates = regmatches(dates, m)
    time(tmp_ras) = strptime(dates, "%Y-%m-%d")
  }
  
  # Y = values(tmp_ras)
  #dim(Y)   = dims[c(2,1,3)]    # Assign column-major dim expected by Rbeast
  #  avoid the need to transpose
  Y = as.matrix(tmp_ras, wide=TRUE)
  dim(Y)   = dims[c(1,2,3)]  #  could just use dims directly...
  
  metadata = list(
    time             = time(tmp_ras),
    isRegularOrdered = FALSE,    # IRREGULAR input
    whichDimIsTime   = 3,        # 437 is the ts length, so set it to '3' here.
    # time$datestr     = datestr,  # date info is contained in the file names
    # time$strfmt      = 'LT05_018032_20080311.yyyy-mm-dd',
    #deltaTime        = 16/365,     # MODIS data are 16 days
    deltaTime        = 1/12,
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
  
  
  o = beast123 (Y, metadata = metadata, extra = extra, ...)
  
  #  pack in some extra info - underhanded but works for now
  o$ext = as.vector(ext(raster)) #  ext objects do not survive serialisation otherwise
  o$crs = crs(raster)
  
  gc()
  return (o)
}

parse_arcgis_coord = function (coord) {
  coord = stri_trim(coord)
  
  in_metres = stri_sub(coord, -1) == "m"
  
  xy = stri_split_regex(coord, "\\s")[[1]][1:2]
  xy = stri_replace_all_fixed (xy, replacement="", pattern=",")
  
  if (!in_metres) {
    library("parzer")
    library("sf")
    x = parse_lon(xy[1])
    y = parse_lat(xy[2])
    point = st_sfc(st_point(c(x,y)), crs = 4326)
    m = r"(PROJCS["modis_sinusoidal",GEOGCS["GCS_Unknown_datum_based_upon_the_custom_spheroid",DATUM["D_Not_specified_based_on_custom_spheroid",SPHEROID["Custom_spheroid",6371007.181,0.0]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",0.0],UNIT["Meter",1.0]])"
    crs_modis_sinusoidal = st_crs(m)
    point2 = st_transform (point, crs_modis_sinusoidal)
    p = st_coordinates(point2)
    return (c(p[1], p[2]))
  }
  
  x_hemi = stri_extract(xy[1], regex = "[NESW]")
  x = stri_replace_all_fixed(xy[1], replacement = "", pattern = x_hemi)
  x = as.numeric(x)
  if (x_hemi == "W") {
    x = x * -1
  }
  
  y_hemi = stri_extract(xy[2], regex = "[NESW]")
  y = stri_replace_all_fixed(xy[2], replacement = "", pattern = y_hemi)
  y = as.numeric(y)
  if (y_hemi == "S") {
    y = y * -1
  }
  
  return (c(x, y))
}


ext_from_arcgis_coord = function (coord, xoff=100000, yoff=-100000) {
  coord = parse_arcgis_coord(coord)
  x1 = coord[1]
  x2 = coord[1] + xoff
  y1 = coord[2]
  y2 = coord[2] + yoff
  ext(min(x1, x2), max(x1, x2), min(y1, y2), max(y1, y2))
}

#  convert a coordinate copied from ArcGIS to a format usable with terra
#  requires the extent be stored on the rbeast object
coord2idx_rbeast = function (b, coord) {
  coord = parse_arcgis_coord(coord)
  x = coord[1]
  y = coord[2]
  
  extent = ext(b$ext)
  
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


plot_beast_modis_coord = function (b, coord, t=FALSE) {
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
  plot (b[rc])
  
}

load_data = function (file = NULL, drivers=NULL) {
  if (is.null(file)) {
    file = file.choose()
  }
  r = rast(file, drivers=drivers)
  if (any (is.na(time(r)))) {
    dates = strptime(names(r), "%Y-%m-%d")
    if (any(is.na(dates))) {
      stop ("Some of the field names do not satisfy the date format requirement (yyyy-mm-dd)")
    }
    time(r) = as.Date(dates)
  }
  return (r)
}


beast_time_to_date = function (x) {
  #  this fails because the current month and day are assumed for %Y on its own
  #b$time_as_date = as.Date(as.character(floor(b$time)), format="%Y") + (b$time - floor (b$time))*365
  as.Date(paste(floor(x), ceiling((x - floor (x))*365), sep=""), format="%Y%j")
}

rasterise_beast = function (b) {
  if (is.null(b$ext)) {
    stop ("An extent must be set on the beast object")
  }
  
  b$time_as_date = lubridate::date_decimal (b$time)
  
  results = list ()
  trend_ncp = beastbit2raster (b, "trend", "ncp")
  
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

beastbit2raster = function (b, component = "trend", subcomponent = "ncp", template=NULL) {
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
        message (sprintf ("Temporal data for %s, %s", component, subcomponent))
        t = b$time_as_date  #  use pre-calculated if exists
        if (is.null(t)) {
          t = lubridate::date_decimal (b$time)
        }
        time(r) = t
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

export_beast_rasters = function (b, dir, prefix="", overwrite=FALSE) {
  list = rasterise_beast(b)
  message ("Exporting now")
  for (name in names(list)) {
    pfx = file.path (dir, paste0(prefix, name))
    outfile = paste0 (pfx, ".tif")
    message (outfile)
    r = list[[name]]
    r[r == Inf] = NA
    writeRaster(list[[name]], outfile, overwrite=overwrite)
    #  also dump netCDF for temporal data
    if (!anyNA(time(r))) {
      outfile = paste0 (pfx, ".nc")
      message (outfile)
      writeCDF(list[[name]], outfile, overwrite=overwrite)
    }
  }
  invisible(list)
}

