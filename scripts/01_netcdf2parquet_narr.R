# Environmental data: process NetCDF files to Parquet for DuckDB analysis

## MUST get through step 2 before running the code!

# 1: Install packages
# repo_line <- 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))'
# writeLines(repo_line, "~/.Rprofile")
# Restart R session to lock in Rprofile changes
# Session -> Restart R or Ctrl+Shift+F10

# For netcdf and raster package (in terminal):
#sudo apt update
#sudo apt install r-cran-ncdf4
#sudo apt-get install libgdal-dev libgeos-dev libproj-dev libsqlite3-dev

# install.packages("terra", repos = "https://rspatial.r-universe.dev")
# install.packages("raster")
# remotes::install_version("raster", version = "3.6-26")

library(ncdf4)
library(raster)
library(lubridate)
library(arrow)
library(dplyr)
library(tidyr)

# 2: Mount the bucket.
# In the terminal, run: bash scripts/00_mount_bucket_folder.sh

# 3: Configuration
# Variables
dataset <- "narr_vwnd"
datatype <- "vwnd" # "air" for air2m and airsfc
years <- c(2025)   # use 2020 or c(2020, 2021, 2022) or c(2004:2025)

# To check variable name for a datatype, if needed
# nc <- ncdf4::nc_open(file)
# names(nc$var)
# print(nc)
# ncdf4::nc_close(nc)

# input_dir <- "~/environ/01_source/narr_acpcp" # mounted bucket location
# output_root <- "~/environ/02_parquet/narr_acpcp" # parquet output root

input_dir <- paste0("P:/environ/01_source/", dataset) # mounted bucket location
output_root <- paste0("P:/environ/02_parquet/", dataset) # parquet output root

# Create output file location
# unlink(output_root, recursive = TRUE) # only use if needing to start over from scratch; otherwise, skip this step since it deletes everything for the variable
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

grid_dir <- file.path(output_root, "grid")
env_dir <- file.path(output_root, "env")

dir.create(grid_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(env_dir, recursive = TRUE, showWarnings = FALSE)

# Target grid / projection used in the old workflow
clip_extent <- extent(-5000000, 1900000, 300000, 4000000)
target_raster <- raster(
  clip_extent,
  resolution = c(25067.53, 25067.53),
  crs = CRS("+init=epsg:3338")
)

# -------------------------
# Build and write grid table once
# -------------------------

grid_xy <- raster::xyFromCell(target_raster, 1:ncell(target_raster))

grid_df <- data.frame(
  cell_id = seq_len(nrow(grid_xy)),
  x = grid_xy[, 1],
  y = grid_xy[, 2],
  crs_epsg = 3338
)

grid_path <- file.path(grid_dir, "grid_3338.parquet")
arrow::write_parquet(grid_df, grid_path)

message("Wrote grid table: ", grid_path)

# -------------------------
# Helper functions
# -------------------------

parse_nc_time <- function(time_vals, time_units) {
  time_origin <- sub(".*since ", "", time_units)
  
  if (grepl("second", time_units, ignore.case = TRUE)) {
    as.POSIXct(time_origin, tz = "UTC") + time_vals
  } else if (grepl("minute", time_units, ignore.case = TRUE)) {
    as.POSIXct(time_origin, tz = "UTC") + time_vals * 60
  } else if (grepl("day", time_units, ignore.case = TRUE)) {
    as.POSIXct(time_origin, tz = "UTC") + time_vals * 86400
  } else {
    as.POSIXct(time_origin, tz = "UTC") + time_vals * 3600
  }
}

raster_to_env_table <- function(r, timestamp, time_index, dataset_name, source_file, cell_ref) {
  pts <- raster::rasterToPoints(r, spatial = FALSE)
  if (is.null(pts) || nrow(pts) == 0) return(NULL)
  
  df <- as.data.frame(pts)
  names(df) <- c("x", "y", "value")
  
  # Match each raster point back to the target grid by coordinates
  # (exact match because r_proj is aligned to target_raster)
  out <- cell_ref |>
    dplyr::left_join(df, by = c("x", "y")) |>
    dplyr::filter(!is.na(value)) |>
    dplyr::mutate(
      dataset = dataset_name,
      timestamp = timestamp,
      time_index = time_index,
      source_file = source_file,
      year = lubridate::year(timestamp),
      month = lubridate::month(timestamp),
      day = lubridate::day(timestamp),
      hour = lubridate::hour(timestamp)
    ) |>
    dplyr::select(dataset, year, month, day, hour, timestamp, time_index, cell_id, x, y, value, source_file)
  
  out
}

# -------------------------
# Process files
# -------------------------

files <- list.files(input_dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)

# keep only files whose names contain one of the requested years
year_pattern <- paste0("(", paste(years, collapse = "|"), ")")
files <- files[grepl(year_pattern, basename(files))]

if (length(files) == 0) stop("No NetCDF files found for requested years in: ", input_dir)

for (file in files) {
  message("Processing file: ", file)

  file_base <- basename(file)

  yr_match <- regmatches(file_base, regexpr("20[0-9]{2}", file_base))
  if (length(yr_match) == 0 || yr_match == "") {
    warning("Could not infer year from filename: ", file_base)
    next
  }
  yr <- as.integer(yr_match)

  if (!(yr %in% years)) next
  
  nc <- ncdf4::nc_open(file)
  
  if (!(datatype %in% names(nc$var))) {
    warning("Variable '", datatype, "' not found in file: ", file_base)
    ncdf4::nc_close(nc)
    next
  }
  
  if (!("time" %in% names(nc$dim))) {
    warning("No time dimension found in file: ", file_base)
    ncdf4::nc_close(nc)
    next
  }
  
  arr <- ncdf4::ncvar_get(nc, datatype)
  dims <- dim(arr)
  
  if (length(dims) != 3) {
    warning("Unexpected dimensions for ", file_base, ": ", paste(dims, collapse = " x "))
    ncdf4::nc_close(nc)
    next
  }
  
  time_vals <- ncdf4::ncvar_get(nc, "time")
  time_units <- ncdf4::ncatt_get(nc, "time", "units")$value
  times <- parse_nc_time(time_vals, time_units)
  
  if (length(times) != dims[3]) {
    warning(
      "Time count does not match slice count for ", file_base,
      " (time=", length(times), ", slices=", dims[3], ")"
    )
  }
  
  chunk_size <- 500
  chunk_rows <- list()
  part_idx <- 1
  
  for (t_idx in seq_len(dims[3])) {
    message("  time slice ", t_idx, " / ", dims[3])
    
    mat <- arr[, , t_idx]
    
    # Build raster from matrix and align to target grid
    r <- raster::raster(t(mat[nrow(mat):1, ]))
    extent(r) <- extent(target_raster)
    crs(r) <- crs(target_raster)
    
    r_proj <- raster::projectRaster(r, target_raster, method = "bilinear")
    
    ts <- if (t_idx <= length(times)) times[t_idx] else as.POSIXct(paste0(yr, "-01-01 00:00:00"), tz = "UTC")
    
    slice_df <- raster_to_env_table(
      r = r_proj,
      timestamp = ts,
      time_index = t_idx,
      dataset_name = dataset,
      source_file = file_base,
      cell_ref = grid_df
    )
    
    if (!is.null(slice_df) && nrow(slice_df) > 0) {
      chunk_rows[[length(chunk_rows) + 1]] <- slice_df
    }
    
    if (length(chunk_rows) >= chunk_size) {
      out_df <- dplyr::bind_rows(chunk_rows)
      
      out_path <- file.path(
        env_dir,
        paste0(tools::file_path_sans_ext(file_base), "_part", sprintf("%03d", part_idx), ".parquet")
      )
      
      arrow::write_parquet(out_df, out_path)
      message("  wrote chunk: ", out_path)
      
      chunk_rows <- list()
      part_idx <- part_idx + 1
    }
  }
  
  # Flush remaining rows
  if (length(chunk_rows) > 0) {
    out_df <- dplyr::bind_rows(chunk_rows)
    
    out_path <- file.path(
      env_dir,
      paste0(tools::file_path_sans_ext(file_base), "_part", sprintf("%03d", part_idx), ".parquet")
    )
    
    arrow::write_parquet(out_df, out_path)
    message("  wrote final chunk: ", out_path)
  }
  
  ncdf4::nc_close(nc)
}