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
    if (i %% 10 == 1) {
      message (sprintf ("%d of %d", i, ncells))
    }
    
    u = unlist(raster[c])

    t2 = bfast::bfastts(u, dates, type = '16-day')
    #decomp = ifelse (na_frac == 0, "stlplus", "stl")
    #  consider suppressWarnings
    tb = bfast::bfast(t2, h=h, decomp="stlplus")
    o[[c]] = tb
  }

  e = ext(raster)
  dims = dim(raster)
  res = list (models = o, ext = e, dims = dims)

  invisible (res)
}


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
