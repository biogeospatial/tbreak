library("terra")
library("sf")
library("stringi")

run_and_save_by_tiles = function (raster, output_dir, output_prefix = "", overwrite_flag=FALSE, tile_size=256, ...) {
  v = get_raster_tile_sf(raster, tile_size=tile_size)

  dir.create(file.path(output_dir), showWarnings = FALSE)

  b = list()
  subset = sample(1:nrow(v))
  #subset = c(111)
  ntiles = nrow(v)

  if (stringi::stri_length(output_prefix) > 0 && !stringi::stri_endswith_fixed(output_prefix, "_") ) {
    output_prefix = paste0(output_prefix, "_")
  }
  #subset = 1:3  #  for debug
  for (i in v$beast_id[subset]) {
    if (is.na(i)) {break}  #  for debug
    pfx = sprintf("%stile%s_", output_prefix, i)
    rds_file = file.path(output_dir, paste0(pfx, ".rds"))

    chk = file.path(output_dir, paste0(pfx, "R2.tif"))
    #message (chk)
    if (file.exists(chk)) {
      next
    }

    message (sprintf("run and save tile %s of %s", i, ntiles))
    r = crop(raster, v[i,], ext=TRUE)

    #  skip if all nodata
    if (is.na(max(minmax(r)))) {
      next
    }

    #devtools::reload("Rbeast")

    beast_result = tbreak:::tiled_beast_modis(r, printParameter=FALSE, printProgress=FALSE, ...)
    #saveRDS(beast_result, rds_file)

    b_rasters = tbreak:::export_beast_rasters(beast_result, output_dir, pfx, overwrite_flag)
    #b[[i]] = b_rasters
    gc()
  }

  invisible(b)
}


get_raster_tile_sf = function (raster, tile_size=512) {
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
  v = sf::st_as_sf(v)  #  make it an sf object

  rm (rr)
  gc()

  return (v)
}
