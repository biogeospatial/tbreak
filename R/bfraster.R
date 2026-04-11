library("bfast")


bfast_raster = function (raster, h=1/7, ...) {

  rdims = dim(raster)

  if (rdims[3] < 20) {
    stop ("Fewer than 20 time steps in data set")
  }


  dates = get_dates_from_raster (raster)

  #  make h proportional to time period
  if (h > 1) {
    h = h / as.numeric (dates[length(dates)] - dates[1])
  }

  targets = (sum(is.na(raster)) / rdims[3]) <= 0.1
  targets[targets == 0] = NA

  #  pre-allocate vector - might save some computation
  o = vector("list", rdims[1] * rdims[2])

  cell_ids = cells(targets)
  ncells = length(cell_ids)

  i = 0
  for (c in cell_ids) {
    i = i + 1
    if (i %% 100 == 1) {
      message (sprintf ("Processing cell %d of %d", i, ncells))
    }

    #  Probably inefficient to get each cell separately but
    #  avoids converting the whole raster to a matrix and
    #  thus doubling the memory load.
    #  Profiling shows the bfast::bfast call takes all the time anyway.
    u = unlist(raster[c])

    t2 = bfast::bfastts(u, dates, type = '16-day')
    #decomp = ifelse (na_frac == 0, "stlplus", "stl")
    #  consider suppressWarnings
    tb = bfast::bfast(t2, h=h, decomp="stlplus")
    o[[c]] = tb
  }

  dims = dim(raster)
  e = as.vector(ext(raster)) #  ext objects do not survive serialisation otherwise
  crs = crs(raster)
  res = list (models = o, ext = e, dims = dims, crs = crs)

  invisible (res)
}

#  very similar to get_date_vector_from_raster
get_dates_from_raster = function (raster){
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
  dates
}

bfast_bit_to_raster = function (b, component="trend_ncp") {

  stopifnot(
    !is.null(component),
    component %in% c("trend_ncp", "season_ncp", "r2")
  )

  funcs = list (
    "trend_ncp" = function (m) {
      ifelse (
        is.null(m),
        NA,
        sum(m$output[[length(m$output)]]$Vt.bp > 0)
      )
    },
    "season_ncp" = function (m) {
      ifelse (
        is.null(m),
        NA,
        sum(m$output[[length(m$output)]]$Wt.bp > 0)
      )
    },
    "r2" = function (m) {
      ifelse (
        is.null(m),
        NA,
        1 - sum (m$output[[length(m$output)]]$Nt ** 2, na.rm=TRUE)
        / sum (m$Yt ** 2, na.rm=TRUE)
      )
    }
  )

  f = funcs[[component]]
  vals = sapply (b$models, FUN = f)

  rr = rast(
    nrows=b$dims[1],
    ncols=b$dims[2],
    crs=b$crs,
    extent=terra::ext(b$ext),
    vals = vals
  )

  return(rr)
}

export_bfast_rasters = function (b, dir, prefix="", overwrite=FALSE) {

  model1 = b$models[[1]]
  stopifnot(class(model1) == "bfast")

  outputs = character()


  for (name in c("trend_ncp", "season_ncp", "r2")) {
    pfx = file.path (dir, paste0(prefix, name))
    outfile = paste0 (pfx, ".tif")
    message (outfile)

    rr = bfast_bit_to_raster(b, name)

    writeRaster(rr, outfile, overwrite=overwrite)
    outputs = c(outputs, outfile)
  }

  invisible(outputs)
}

