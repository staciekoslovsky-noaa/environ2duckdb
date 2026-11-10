# To install RPostgres, in terminal
# sudo apt-get install libpq-dev libssl-dev

# Load libraries
library(RPostgres)
library(tidyverse)
library(sf)
library(raster)
library(duckdb)

# Set variables -----------------------------------------
# Set start/end date 
start_date <- '2025-04-04'
end_date <- '2025-06-11'

# Select NARR variables
narr_variables <- list(
  "narr_acpcp",
  "narr_air2m",
  "narr_airsfc",
  "narr_prmsl",
  "narr_uwnd",
  "narr_vwnd"
)

# Indicate if sea ice data are provided
sea_ice <- TRUE # indicates whether to include an export for "nsidc_cdr_v6" data

# Process the grid from the database to duckDB
con <- RPostgres::dbConnect(RPostgres::Postgres(),
                              dbname = Sys.getenv("mml_pep_db"), 
                              host = Sys.getenv("mml_pep_ip"), 
                              user = Sys.getenv("mml_pep_user"), 
                              password = Sys.getenv("mml_pep_user_pw"))

centroids <- sf::st_read(con, layer = Id(schema = "base", table  = "geo_analysis_grid_centroid"))

target_resolution <- 25067.53
grid <- terra::rast(
  ncols = round((1900000 - (-5000000)) / target_resolution),
  nrows = round((4000000 - 300000) / target_resolution),
  xmin = -5000000,
  xmax = 1900000,
  ymin = 300000,
  ymax = 4000000,
  crs = "EPSG:3338"
)

centroids_xy <- sf::st_coordinates(sf::st_transform(centroids, 3338))

centroids <- sf::st_drop_geometry(centroids)
centroids$cell_id <- raster::cellFromXY(grid, centroids_xy[, c("X", "Y")])
centroids$input_row <- seq_len(nrow(centroids))

duck <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
DBI::dbExecute(duck, "INSTALL spatial;")
DBI::dbExecute(duck, "LOAD spatial;")

duckdb::duckdb_register(duck, "centroids", centroids, overwrite = TRUE)
rm(grid, centroids_xy)

# Process selected covariates -----------------------------------------
get_covariates <- function(variable, start_date, end_date, duck) {
  parquet_path <- paste0("P:\\environ\\02_parquet\\", variable, "\\env") # parquet output root
  path <- DBI::dbQuoteString(duck, path.expand(parquet_path))

  sql <- paste0(
  "SELECT c.cell, e.timestamp AS covariate_dt, e.value AS ", variable,
  " FROM centroids c ",
  "LEFT JOIN read_parquet(", path, ", hive_partitioning = false) e ",
  "ON c.cell_id = e.cell_id ",
  "WHERE e.timestamp >= \'", start_date,
  "\'::TIMESTAMPTZ AND e.timestamp <= \'", end_date,
  "\'::TIMESTAMPTZ"
  )

  covariates <- DBI::dbGetQuery(duck, sql) %>%
    mutate(cell = as.integer(cell)) 
}

# Process NARR data
for (i in 1:length(narr_variables)) {
  narr_variable <- narr_variables[i]
  data_i <- get_covariates(narr_variable, start_date, end_date, duck)

  if (i == 1) {
    covariates_narr <- centroids %>%
      inner_join(data_i, by = "cell")
    rm(data_i)
  } else {
    covariates_narr <- covariates_narr %>%
      inner_join(data_i, by = c("cell", "covariate_dt"))
    rm(data_i)
  }

}
covariates_narr <- covariates_narr %>%
  select(-cell_id, -input_row)

# Process sea ice data (not available for all grid cells due to land definition)
covariates_nsidc <- get_covariates("nsidc_cdr_v6", start_date, end_date, duck)

# Clean up workspace
DBI::dbDisconnect(duck)
rm(i, duck, narr_variable)

# Write data to DB
target_narr <- Id(
    schema = "surv_ice_seals_2025",
    table = "tbl_env_narr"
  )

 RPostgres::dbWriteTable(
    con,
    target_narr,
    covariates_narr,
    append = FALSE,
    row.names = FALSE
  )

target_nsidc <- Id(
    schema = "surv_ice_seals_2025",
    table = "tbl_env_nsidc"
  )

  RPostgres::dbWriteTable(
    con,
    target_nsidc,
    covariates_nsidc,
    append = FALSE,
    row.names = FALSE
  )

RPostgres::dbDisconnect(con)
rm(con)


# SELECT *
# FROM read_parquet(
#   'P:/environ/02_parquet/nsidc_cdr_v6/env/nsidc_cdr_v6_*.parquet'
# );