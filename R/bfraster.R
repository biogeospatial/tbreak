library("bfast")
library("foreach")

plot_bfast_modis_coord = function (raster, coord, h=0.15, main=NULL) {
  coord = parse_coord_string(coord, crs(raster))

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
  #u[u < -0.25] = NA

  na_frac = sum(is.na(u)) / length(u)
  if (na_frac > 0.2) {
    message (
      sprintf (
        "More than 20% of records are NA (%d), skipping bfast generation",
        na_frac * 100
      )
    )
    return()
  }

  #  proportional to time period
  if (h > 1) {
    h = h / as.numeric (dates[length(dates)] - dates[1])
  }

  t2 = bfast::bfastts(u, dates, type = '16-day')
  tb = bfast::bfast(t2, h=h)
  plot(tb, main=main)
  invisible (tb)
}

bfast_raster = function (raster, h=1/7, dopar=FALSE, ...) {

  rdims = dim(raster)

  if (rdims[3] < 20) {
    stop ("Fewer than 20 time steps in data set")
  }


  dates = get_dates_from_raster (raster)

  #  make h proportional to time period
  if (h > 1) {
    h = h / as.numeric (dates[length(dates)] - dates[1])
  }

  targets = (sum(is.na(raster)) / rdims[3]) <= 0.2
  targets[targets == 0] = NA

  #  pre-allocate vector - might save some computation
  o = vector("list", rdims[1] * rdims[2])

  cell_ids = cells(targets)
  ncells = length(cell_ids)

  if (!dopar) {
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

      tryCatch ({
        t2 = bfast::bfastts(u, dates, type = '16-day')
        #decomp = ifelse (na_frac == 0, "stlplus", "stl")
        #  consider suppressWarnings
        tb = bfast::bfast(t2, h=h, decomp="stlplus")
      },
      warning = function (w) {
        w_coord = xyFromCell(raster, c)
        warning (sprintf ("Coord %s %s: %s", w_coord[1], w_coord[2], w))
      },
      error = function (e) {
        tb = NULL
      })

      o[[c]] = tb
    }
  } else {

    library(doParallel)
    cl <- makeCluster(detectCores() - 1) # Leaves one core free for OS tasks
    registerDoParallel(cl)

    #  don't enclose raster objects in the function or we get nasty errors
    #  prob due to gc.
    cells = lapply (cell_ids, FUN = function (c) { raster[c] })

    oo = foreach (c=cells, .packages=c("bfast")) %dopar% {
      u = unlist(c)
      tryCatch ({
        t2 = bfast::bfastts(u, dates, type = '16-day')
        tb = bfast::bfast(t2, h=h, decomp="stlplus")
      },
      warning = function (e) {
        warning (sprintf ("Cell %s: %s", c, e))
      },
      error = function (e) {
        tb = NULL
      })
    }

    stopCluster(cl)


    #  now populate the main list
    if (length(cell_ids) != rdims[1] * rdims[2]) {
      for (i in 1:length(cell_ids)) {
        o[[cell_ids[i]]] = oo[[i]]
      }
    } else {
      o = oo
    }
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
    component %in% c("trend_ncp", "season_ncp", "R2")
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
    "R2" = function (m) {
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


  for (name in c("trend_ncp", "season_ncp", "R2")) {
    pfx = file.path (dir, paste0(prefix, name))
    outfile = paste0 (pfx, ".tif")
    message (outfile)

    rr = bfast_bit_to_raster(b, name)

    writeRaster(rr, outfile, overwrite=overwrite)
    outputs = c(outputs, outfile)
  }

  invisible(outputs)
}

