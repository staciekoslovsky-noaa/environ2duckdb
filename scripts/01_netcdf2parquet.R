# Environ_01_LoadNARR_acpcp_parquet.R
# Convert NARR ACPCP NetCDF files to Parquet for DuckDB analysis

# Before running this code, need to mount the bucket.
# In the terminal, run: bash scripts/00_mount_bucket_folder.sh

install_pkg <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    if (!require(x, character.only = TRUE)) stop("Package not found: ", x)
  }
}

install_pkg("ncdf4")
install_pkg("raster")
install_pkg("lubridate")
install_pkg("arrow")
install_pkg("dplyr")
install_pkg("tidyr")

library(ncdf4)
library(raster)
library(lubridate)
library(arrow)
library(dplyr)
library(tidyr)

# ---- configuration ----
input_dir <- "/home/staciekoslovsky/my_gcs_bucket/01_source/narr_acpcp" # mounted bucket location
output_dir <- "/home/staciekoslovsky/my_gcs_bucket/02_parquet" # parquet output root
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# target grid used in the old workflow
clip_extent <- extent(-5000000, 1900000, 300000, 4000000)
target_raster <- raster(
  clip_extent,
  resolution = c(25067.53, 25067.53),
  crs = CRS("+init=epsg:3338")
)

# helper: convert raster to tidy table
raster_to_table <- function(r, timestamp, source_file, band_index) {
  pts <- rasterToPoints(r, spatial = FALSE)
  if (is.null(pts) || nrow(pts) == 0) return(NULL)

  df <- as.data.frame(pts)
  names(df) <- c("x", "y", "value")

  df %>%
    filter(!is.na(value)) %>%
    mutate(
      timestamp = timestamp,
      source_file = source_file,
      band = band_index,
      year = year(timestamp),
      month = month(timestamp),
      day = day(timestamp),
      hour = hour(timestamp),
      crs_epsg = 3338
    ) %>%
    select(year, month, day, hour, timestamp, band, x, y, value, crs_epsg, source_file)
}

# list NetCDF files
files <- list.files(input_dir, pattern = "\\.nc$", full.names = TRUE)

# process each file
for (file in files) {
  message("Processing: ", file)

  # infer year from filename if needed
  # adjust this if filenames use another convention
  file_base <- basename(file)
  year_match <- regmatches(file_base, regexpr("20[0-9]{2}", file_base))
  if (length(year_match) == 0 || year_match == "") {
    warning("Could not infer year from filename: ", file_base)
    next
  }
  yr <- as.integer(year_match)

  nc <- nc_open(file)
  on.exit(nc_close(nc), add = TRUE)

  # identify main variable if needed
  var_names <- names(nc$var)
  if (length(var_names) == 0) {
    warning("No variables found in file: ", file_base)
    nc_close(nc)
    next
  }

  var_name <- var_names[1]  # adjust if ACPCP variable is named explicitly
  var_data <- ncvar_get(nc, var_name)

  # if file has time dimension, pull time values
  time_var <- NULL
  if ("time" %in% names(nc$dim)) {
    time_vals <- ncvar_get(nc, "time")
    time_units <- ncatt_get(nc, "time", "units")$value
    time_origin <- sub(".*since ", "", time_units)
    time_var <- as.POSIXct(time_origin, tz = "UTC") + time_vals * 3600
  }

  # determine band count
  if (length(dim(var_data)) == 3) {
    n_bands <- dim(var_data)[3]
  } else if (length(dim(var_data)) == 2) {
    n_bands <- 1
  } else {
    warning("Unexpected variable dimensions in file: ", file_base)
    nc_close(nc)
    next
  }

  all_rows <- list()

  for (b in seq_len(n_bands)) {
    message("  band ", b, " of ", n_bands)

    if (n_bands == 1) {
      band_mat <- var_data
    } else {
      band_mat <- var_data[, , b]
    }

    # convert matrix to raster
    r <- raster(t(band_mat[nrow(band_mat):1, ]))
    extent(r) <- extent(target_raster)
    crs(r) <- crs(target_raster)

    # project to target grid
    r_proj <- projectRaster(r, target_raster, method = "bilinear")

    ts <- if (!is.null(time_var) && length(time_var) >= b) time_var[b] else as.POSIXct(paste0(yr, "-01-01 00:00:00"), tz = "UTC")

    band_df <- raster_to_table(
      r_proj,
      timestamp = ts,
      source_file = file_base,
      band_index = b
    )

    if (!is.null(band_df)) {
      all_rows[[length(all_rows) + 1]] <- band_df
    }
  }

  if (length(all_rows) > 0) {
    out_df <- bind_rows(all_rows)

    out_path <- file.path(output_dir, paste0("acpcp_", yr, "_", tools::file_path_sans_ext(file_base), ".parquet"))
    write_parquet(out_df, out_path)
    message("Wrote: ", out_path)
  }

  nc_close(nc)
}