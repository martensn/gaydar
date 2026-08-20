# code/build_national_tract_rates.R
#
# One-time offline build: combines the 51 already-cached per-state tract
# layers (data/cache/tract_state/*.rds) into a single lightweight,
# geometry-free national lookup table. Used at runtime for state/national
# percentile comparisons (app.R) and as the source of each tract's
# Congressional District assignment for the geography selector.
#
# No population re-estimation happens here -- this only reads the existing
# tract-level rates and joins in a Congressional District id. The only
# network calls are fetching CD cartographic boundaries via tigris (cached
# locally after first use, same as the rest of the pipeline).
#
# Run once locally: Rscript code/build_national_tract_rates.R

library(dplyr)
library(sf)
library(purrr)
library(stringr)

# Planar GEOS instead of spherical S2 for the CD overlap join below: S2 is
# stricter about topology and can error on a degenerate edge (duplicate
# vertex) that GEOS tolerates fine. This script only compares relative
# overlap areas within a state to pick the largest match, so the small
# planar-vs-geodesic difference doesn't matter here (unlike the antimeridian
# handling in code/helpers.R's IDW step, which genuinely needs S2).
sf::sf_use_s2(FALSE)

root_dir  <- "."
cache_dir <- file.path(root_dir, "data/cache/tract_state")
out_path  <- file.path(root_dir, "data/cache/national_tract_rates.rds")

files <- list.files(cache_dir, pattern = "^tract_.*\\.rds$", full.names = TRUE)
stopifnot(length(files) > 0)

state_from_file <- function(f) {
  stringr::str_match(basename(f), "^tract_([A-Z]{2})_")[, 2]
}

keep_cols <- c(
  "GEOID", "puma_id", "total_pop", "total_m", "total_w", "total_nb",
  "lgbt_m_map", "lgbt_w_map", "lgbt_nb_map",
  "lg_m_map",   "lg_w_map",   "lg_nb_map",
  "bi_m_map",   "bi_w_map",   "bi_nb_map",
  "queer_m_map","queer_w_map","queer_nb_map",
  "trans_m_map","trans_w_map"
)

# --- National Congressional District boundaries (single download, once) ---
# tigris::congressional_districts() hardcodes a year -> congress-number
# lookup table that doesn't cover 2023+ (its GENZ2022 entry points at the
# 116th congress, which 404s -- GENZ2022 cb files are actually keyed to the
# 118th). Rather than depend on that stale mapping, fetch the correct file
# directly; it's a single national file (not per-state), so this is one
# download for all 51 states rather than 51.
cd_cache_dir <- file.path(root_dir, "data/cache/cd_boundaries")
dir.create(cd_cache_dir, recursive = TRUE, showWarnings = FALSE)
cd_zip_path <- file.path(cd_cache_dir, "cb_2022_us_cd118_500k.zip")
cd_zip_url  <- "https://www2.census.gov/geo/tiger/GENZ2022/shp/cb_2022_us_cd118_500k.zip"

if (!file.exists(cd_zip_path)) {
  message("Downloading national Congressional District boundaries...")
  download.file(cd_zip_url, cd_zip_path, mode = "wb", quiet = TRUE)
}

cds_national <- sf::st_read(paste0("/vsizip/", cd_zip_path), quiet = TRUE) |>
  dplyr::transmute(cd_id = GEOID, state_fips = STATEFP) |>
  sf::st_make_valid() |>
  sf::st_transform(4326)

# Assigns each tract the Congressional District with greatest overlap area.
# Falls back to all-NA cd_id for a state if the join itself errors, so one
# bad geometry doesn't abort the whole 51-state batch.
assign_cd <- function(tract_sf, state_fips) {
  cds_state <- cds_national |> dplyr::filter(state_fips == !!state_fips)

  if (nrow(cds_state) == 0) {
    return(tibble::tibble(GEOID = tract_sf$GEOID, cd_id = NA_character_))
  }

  tryCatch(
    sf::st_join(
      tract_sf |> dplyr::select(GEOID) |> sf::st_make_valid(),
      cds_state |> dplyr::select(cd_id),
      join = sf::st_intersects,
      largest = TRUE
    ) |>
      sf::st_drop_geometry() |>
      dplyr::distinct(GEOID, .keep_all = TRUE),
    error = function(e) {
      message("    CD join failed (", conditionMessage(e), "), leaving cd_id NA")
      tibble::tibble(GEOID = tract_sf$GEOID, cd_id = NA_character_)
    }
  )
}

message("Building national tract rates table from ", length(files), " state caches...")

national <- purrr::map_dfr(files, function(f) {
  state_abbr <- state_from_file(f)
  message("  ", state_abbr)

  tract_sf <- readRDS(f)
  tract_sf <- tract_sf[sf::st_geometry_type(tract_sf) %in% c("POLYGON", "MULTIPOLYGON"), ]
  tract_sf <- sf::st_make_valid(tract_sf)
  tract_sf <- tract_sf[sf::st_is_valid(tract_sf) & !sf::st_is_empty(tract_sf), ]

  state_fips <- substr(tract_sf$GEOID[1], 1, 2)
  cd_map <- assign_cd(tract_sf, state_fips)

  tract_sf |>
    sf::st_drop_geometry() |>
    dplyr::select(dplyr::any_of(keep_cols)) |>
    dplyr::left_join(cd_map, by = "GEOID") |>
    dplyr::mutate(
      state       = state_abbr,
      county_fips = substr(GEOID, 1, 5)
    )
})

n_missing_cd <- sum(is.na(national$cd_id))
message("Rows: ", nrow(national), " | States: ", dplyr::n_distinct(national$state),
        " | Missing cd_id: ", n_missing_cd, " (", round(100 * n_missing_cd / nrow(national), 2), "%)")

saveRDS(national, out_path)
message("Saved: ", out_path, " (", round(file.size(out_path) / 1e6, 2), " MB)")
