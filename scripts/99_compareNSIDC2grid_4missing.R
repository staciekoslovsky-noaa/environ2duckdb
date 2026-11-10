file <- "P:/environ/01_source/nsidc_cdr_v6/2025/sic_psn25_20250416_am2_v06r00.nc"

# Check the actual variable metadata, including scale and flag definitions.
nc <- ncdf4::nc_open(file)
print(nc$var[["cdr_seaice_conc"]])
print(ncdf4::ncatt_get(nc, "cdr_seaice_conc"))
ncdf4::nc_close(nc)

source <- terra::rast(file, subds = "cdr_seaice_conc")[[1]]

print(source)
print(terra::global(source, range, na.rm = TRUE))

target_raster <- terra::rast(
  ncols = round((1900000 - (-5000000)) / target_resolution),
  nrows = round((4000000 - 300000) / target_resolution),
  xmin = -5000000,
  xmax = 1900000,
  ymin = 300000,
  ymax = 4000000,
  crs = "EPSG:3338"
)

clean <- terra::ifel(source >= 0 & source <= 1, source, NA)
projected <- terra::project(
  clean,
  target_raster,
  method = "bilinear"
)

names(projected) <- "value"

# Unlike the production output, retain every target cell.
diagnostic <- terra::as.data.frame(
  projected,
  cells = TRUE,
  xy = TRUE,
  na.rm = FALSE
)

table(missing = is.na(diagnostic$value))

terra::plot(projected, main = "Projected NSIDC concentration")
terra::plot(is.na(projected), main = "Missing target cells")








con <- RPostgres::dbConnect(RPostgres::Postgres(),
                              dbname = Sys.getenv("mml_pep_db"), 
                              host = Sys.getenv("mml_pep_ip"), 
                              user = Sys.getenv("mml_pep_user"), 
                              password = Sys.getenv("mml_pep_user_pw"))

centroids <- sf::st_read(con, layer = Id(schema = "base", table  = "geo_analysis_grid_centroid"))
RPostgres::dbDisconnect(con)


# Use the single-date projected raster created earlier.
sea_ice_layer <- projected[[1]]

# Transform points into the raster's CRS.
points_projected <- sf::st_transform(
  centroids, 
  terra::crs(sea_ice_layer)
)

# Extract the value of the raster cell containing each centroid.
values <- terra::extract(
  sea_ice_layer,
  terra::vect(points_projected),
  method = "simple",
  ID = FALSE
)

comparison <- points_projected
comparison$sea_ice <- values[[1]]
comparison$status <- ifelse(
  is.na(comparison$sea_ice),
  "Missing",
  "Available"
)

table(comparison$status)

# Overlay centroids on the projected concentration layer.
terra::plot(sea_ice_layer, main = "NSIDC concentration and study centroids")

terra::points(
  terra::vect(comparison),
  pch = 21,
  bg = ifelse(comparison$status == "Missing", "red", "dodgerblue"),
  col = "black",
  cex = 0.6
)

legend(
  "bottomleft",
  legend = c("Available", "Missing"),
  pt.bg = c("dodgerblue", "red"),
  col = "black",
  pch = 21,
  bty = "n"
)

# Keep missing centroids for further inspection.
missing_centroids <- comparison[
  comparison$status == "Missing",
]