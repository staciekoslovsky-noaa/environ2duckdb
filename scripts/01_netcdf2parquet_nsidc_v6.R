# ============================================================
# NSIDC Sea Ice Concentration CDR V6: NetCDF -> annual Parquet
#
# Install these packages once, if needed:
# install.packages(c("ncdf4", "terra", "lubridate", "arrow"))
#
# Mount your bucket before running, if applicable:
# bash scripts/00_mount_bucket_folder.sh
#
# Expected input:
# input_dir/2020/date.nc
# input_dir/2021/date.nc
#
# Outputs:
# output_root/grid/grid_3338.parquet
# output_root/env/nsidc_cdr_v6_2020.parquet
# output_root/env/nsidc_cdr_v6_2021.parquet
# ============================================================

library(ncdf4)
library(terra)
library(lubridate)
library(arrow)

# -------------------------
# 1. Configuration
# -------------------------

dataset <- "nsidc_cdr_v6"
datatype <- "cdr_seaice_conc"
years <- 2004:2025

# Use Northern Hemisphere daily files for the Alaska target grid.
input_dir <- paste0("P:/environ/01_source/", dataset) # mounted bucket location
output_root <- paste0("P:/environ/02_parquet/", dataset) # parquet output root

grid_dir <- file.path(output_root, "grid")
env_dir <- file.path(output_root, "env")

year_dirs <- file.path(input_dir, as.character(years))
missing_dirs <- year_dirs[!dir.exists(year_dirs)]

if (length(missing_dirs) > 0L) {
  stop(
    "Missing input directories; check the mount and requested years: ",
    paste(missing_dirs, collapse = ", ")
  )
}

for (directory in c(grid_dir, env_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(directory)) {
    stop("Could not create output directory: ", directory)
  }
}

# Existing annual outputs are protected from accidental replacement.
annual_paths <- file.path(
  env_dir,
  paste0(dataset, "_", years, ".parquet")
)

if (any(file.exists(annual_paths))) {
  stop(
    "Annual outputs already exist. Move or remove the files you intend ",
    "to regenerate before running: ",
    paste(annual_paths[file.exists(annual_paths)], collapse = ", ")
  )
}

# -------------------------
# 2. Target grid
# -------------------------

# Preserve the original extent with approximately 25-km cells.
# The exact cell size is adjusted to fit this extent.
target_resolution <- 25067.53

target_raster <- terra::rast(
  ncols = round((1900000 - (-5000000)) / target_resolution),
  nrows = round((4000000 - 300000) / target_resolution),
  xmin = -5000000,
  xmax = 1900000,
  ymin = 300000,
  ymax = 4000000,
  crs = "EPSG:3338"
)

grid_xy <- terra::xyFromCell(
  target_raster,
  seq_len(terra::ncell(target_raster))
)

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
# 3. Time helpers
# -------------------------

parse_nc_time <- function(time_vals, time_units) {
  parts <- regmatches(
    time_units,
    regexec(
      "^\\s*(seconds?|minutes?|hours?|days?)\\s+since\\s+(.+)\\s*$",
      time_units,
      ignore.case = TRUE
    )
  )[[1]]

  if (length(parts) != 3L) {
    stop("Unsupported NetCDF time units: ", time_units)
  }

  unit <- sub("s$", "", tolower(parts[2]))

  origin <- lubridate::ymd_hms(
    trimws(parts[3]),
    tz = "UTC",
    truncated = 3,
    quiet = TRUE
  )

  if (is.na(origin)) {
    stop("Could not parse NetCDF time origin: ", parts[3])
  }

  multiplier <- switch(
    unit,
    second = 1,
    minute = 60,
    hour = 3600,
    day = 86400
  )

  times <- origin + as.numeric(time_vals) * multiplier

  if (length(times) == 0L || any(!is.finite(as.numeric(times)))) {
    stop("Missing or invalid NetCDF time values.")
  }

  times
}

read_nc_times <- function(file, variable) {
  nc <- ncdf4::nc_open(file)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  if (!(variable %in% names(nc$var))) {
    stop("Variable '", variable, "' not found in: ", file)
  }

  variable_dims <- vapply(
    nc$var[[variable]]$dim,
    function(x) x$name,
    character(1)
  )

  if (!("time" %in% variable_dims)) {
    stop("Variable has no time dimension in: ", file)
  }

  units <- ncdf4::ncatt_get(nc, "time", "units")
  calendar <- ncdf4::ncatt_get(nc, "time", "calendar")

  if (!isTRUE(units$hasatt)) {
    stop("Missing time units in: ", file)
  }

  if (
    isTRUE(calendar$hasatt) &&
    !(tolower(calendar$value) %in%
      c("standard", "gregorian", "proleptic_gregorian"))
  ) {
    stop("Unsupported calendar '", calendar$value, "' in: ", file)
  }

  parse_nc_time(
    ncdf4::ncvar_get(nc, "time"),
    units$value
  )
}

# -------------------------
# 4. Annual processing
# -------------------------

process_year <- function(yr) {
  year_dir <- file.path(input_dir, as.character(yr))

  files <- sort(list.files(
    year_dir,
    pattern = "\\.nc$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  ))

  if (length(files) == 0L) {
    stop("No NetCDF files found in: ", year_dir)
  }

  out_path <- file.path(
    env_dir,
    paste0(dataset, "_", yr, ".parquet")
  )

  if (file.exists(out_path)) {
    stop("Output already exists: ", out_path)
  }

  # Write beside the final output, then rename only after success.
  temp_path <- tempfile(
    pattern = paste0(dataset, "_", yr, "_"),
    tmpdir = env_dir,
    fileext = ".parquet"
  )

  writer <- NULL
  sink <- NULL
  total_rows <- 0
  seen_times <- numeric()

  on.exit({
    tryCatch(
      {
        if (!is.null(writer)) writer$Close()
      },
      finally = {
        tryCatch(
          {
            if (!is.null(sink)) sink$close()
          },
          finally = {
            if (file.exists(temp_path)) {
              if (unlink(temp_path) != 0L) {
                warning("Could not remove temporary file: ", temp_path)
              }
            }
          }
        )
      }
    )
  }, add = TRUE)

  message("Processing year: ", yr, " (", length(files), " files)")

  for (file in files) {
    message("  Reading: ", basename(file))

    file_base <- basename(file)
    times <- read_nc_times(file, datatype)

    if (any(lubridate::year(times) != yr)) {
      stop("NetCDF timestamps do not match year directory: ", file)
    }

    numeric_times <- as.numeric(times)

    # Avoid double-counting overlapping files or mixed hemispheres.
    if (
      anyDuplicated(numeric_times) > 0L ||
      any(numeric_times %in% seen_times)
    ) {
      stop(
        "Duplicate timestamps detected in: ", file,
        ". Include only one daily product and hemisphere."
      )
    }

    seen_times <- c(seen_times, numeric_times)

    # Read the native NetCDF georeferencing; do not assign EPSG:3338 here.
    source_raster <- terra::rast(file, subds = datatype)

    if (!nzchar(terra::crs(source_raster))) {
      stop("No source CRS could be read from: ", file)
    }

    if (terra::nlyr(source_raster) != length(times)) {
      stop(
        "Raster layer/time count mismatch in ", file,
        ": layers=", terra::nlyr(source_raster),
        ", times=", length(times)
      )
    }

    for (t_idx in seq_len(terra::nlyr(source_raster))) {
      timestamp <- times[t_idx]
      r <- source_raster[[t_idx]]

      # terra decodes the NetCDF scale/offset automatically.
      # V6 concentration is a fraction: 0 = ice-free, 1 = 100% ice.
      r <- terra::ifel(r >= 0 & r <= 1, r, NA)

      r_proj <- terra::project(
        r,
        target_raster,
        method = "bilinear"
      )

      names(r_proj) <- "value"

      pts <- terra::as.data.frame(
        r_proj,
        cells = TRUE,
        xy = FALSE,
        na.rm = TRUE
      )

      n <- nrow(pts)

      if (n == 0L) {
        warning(
          "No valid concentrations on the target grid: ",
          file_base, ", time slice ", t_idx
        )
        next
      }

      cells <- grid_df[pts$cell, , drop = FALSE]

      out_df <- data.frame(
        dataset = rep(dataset, n),
        year = rep(as.integer(lubridate::year(timestamp)), n),
        month = rep(as.integer(lubridate::month(timestamp)), n),
        day = rep(as.integer(lubridate::day(timestamp)), n),
        hour = rep(as.integer(lubridate::hour(timestamp)), n),
        timestamp = rep(timestamp, n),
        time_index = rep(as.integer(t_idx), n),
        cell_id = cells$cell_id,
        x = cells$x,
        y = cells$y,
        value = pts$value,
        source_file = rep(file_base, n)
      )

      batch <- arrow::Table$create(out_df)

      if (is.null(writer)) {
        sink <- arrow::FileOutputStream$create(temp_path)

writer <- arrow::ParquetFileWriter$create(
  schema = batch$schema,
  sink = sink,
  properties = arrow::ParquetWriterProperties$create(
    column_names = names(out_df),
    compression = "snappy"
  )
)
      }

      # Add this slice as a row group in the same annual file.
      writer$WriteTable(batch, chunk_size = n)
      total_rows <- total_rows + n
    }
  }

  if (is.null(writer)) {
    stop("No valid concentrations found for year: ", yr)
  }

  writer$Close()
  writer <- NULL

  sink$close()
  sink <- NULL

  if (!file.rename(temp_path, out_path)) {
    stop("Could not publish annual Parquet file: ", out_path)
  }

  message("Wrote annual file: ", out_path, " (", total_rows, " rows)")
}

# -------------------------
# 5. Run
# -------------------------

for (yr in years) {
  process_year(yr)
}

message("Finished processing all requested years.")