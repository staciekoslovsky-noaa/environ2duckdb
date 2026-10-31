# To install RPostgres, in terminal
# sudo apt-get install libpq-dev libssl-dev

# Load libraries
library(RPostgres)
library(sf)
library(raster)
library(duckdb)

# Prepare for running code
parquet_path <- "P:\\environ\\02_parquet\\dataset=acpcp" # parquet output root

con <- RPostgres::dbConnect(RPostgres::Postgres(),
                              dbname = Sys.getenv("mml_pep_db"), 
                              host = Sys.getenv("mml_pep_ip"), 
                              user = Sys.getenv("mml_pep_user"), 
                              password = Sys.getenv("mml_pep_user_pw"))

centroids <- sf::st_read(con, layer = Id(schema = "base", table  = "geo_analysis_grid_centroid"))
grid <- raster::raster(
  raster::extent(-5000000, 1900000, 300000, 4000000),
  resolution = c(25067.53, 25067.53),
  crs = CRS("epsg:3338")
) # as defined by environmental covariate code

centroids_xy <- sf::st_coordinates(sf::st_transform(centroids, 3338))

centroids <- sf::st_drop_geometry(centroids)
centroids$cell_id <- raster::cellFromXY(grid, centroids_xy[, c("X", "Y")])
centroids$input_row <- seq_len(nrow(centroids))


duck <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
DBI::dbExecute(duck, "INSTALL spatial;")
DBI::dbExecute(duck, "LOAD spatial;")

duckdb::duckdb_register(duck, "centroids", centroids, overwrite = TRUE)

path <- DBI::dbQuoteString(duck, path.expand(parquet_path))

sql <- paste0(
  "SELECT c.*, e.dataset, e.timestamp AS covariate_time, e.value ",
  "FROM centroids c ",
  "LEFT JOIN read_parquet(", path, ", hive_partitioning = false) e ",
  "ON c.cell_id = e.cell_id ",
  "WHERE month >= 3 AND month <= 6 AND year = 2025"
)

covariates <- DBI::dbGetQuery(duck, sql) 
#took less than 10 seconds to get acpcp value (based on the criteria above!)
