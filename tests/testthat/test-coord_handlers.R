source("../../R/coord_handlers.R")

describe("parse_coord_string", {

  it("parses nearmap lat/lon", {
    skip_if_not_installed("parzer")
    result <- parse_coord_string("-33, 150")
    expect_length(result, 2)
    expect_equal(result[1], 150)
    expect_equal(result[2], -33)
  })

  it("parses QGIS lat/lon", {
    skip_if_not_installed("parzer")
    result <- parse_coord_string("-33,150")
    expect_length(result, 2)
    expect_equal(result[1], 150)
    expect_equal(result[2], -33)
  })

  it("parses ArcGIS lat/lon with degree symbols", {
    skip_if_not_installed("parzer")
    result <- parse_coord_string("150.7572609°E 33.9385950°S ")
    expect_length(result, 2)
    expect_equal(result[1], 150.7572609, tolerance = 1e-5)
    expect_equal(result[2], -33.9385950, tolerance = 1e-5)
  })

  it("parses ArcGIS projected coordinates with m suffix", {
    result <- parse_coord_string("13,907,112.15 -3,773,782.94 m ")
    expect_length(result, 2)
    expect_equal(result[1], 13907112.15)
    expect_equal(result[2], -3773782.94)
  })

  it("parses QGIS projected coordinates", {
    result <- parse_coord_string("13907112.15,-3773782.94")
    expect_length(result, 2)
    expect_equal(result[1], 13907112.15)
    expect_equal(result[2], -3773782.94)
  })

  it("parses Google Earth DMS format", {
    skip_if_not_installed("parzer")
    result <- parse_coord_string("33°56'19.03\"S 150°45'18.22\"E")
    expect_length(result, 2)
    expect_equal(result[1], 150.75506, tolerance = 1e-4)
    expect_equal(result[2], -33.93862, tolerance = 1e-4)
  })

})

describe("coord2idx_rbeast", {

  it("returns correct row and column for a projected coordinate", {
    skip_if_not_installed("terra")
    b <- list(
      crs   = NULL,
      ext   = c(13900000, 13910000, -3780000, -3770000),
      ncols = 10,
      nrows = 10
    )
    result <- coord2idx_rbeast(b, "13907112.15,-3773782.94")
    expect_length(result, 2)
    expect_equal(result[1], 4)  # row
    expect_equal(result[2], 8)  # col
  })

  it("errors when x coordinate is outside extent", {
    skip_if_not_installed("terra")
    b <- list(
      crs   = NULL,
      ext   = c(13900000, 13910000, -3780000, -3770000),
      ncols = 10,
      nrows = 10
    )
    expect_error(
      coord2idx_rbeast(b, "13800000,-3773782.94"),
      "X coord is outside data set extent"
    )
  })

  it("errors when y coordinate is outside extent", {
    skip_if_not_installed("terra")
    b <- list(
      crs   = NULL,
      ext   = c(13900000, 13910000, -3780000, -3770000),
      ncols = 10,
      nrows = 10
    )
    expect_error(
      coord2idx_rbeast(b, "13907112.15,-3900000"),
      "Y coord is outside data set extent"
    )
  })

})
