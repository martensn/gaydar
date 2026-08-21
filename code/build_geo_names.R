# code/build_geo_names.R
#
# One-time offline build: human-readable names for county/PUMA/Congressional
# District GEOIDs, so the map's tooltip/popup can show e.g. "Alger County, MI"
# instead of "26003". Kept separate from build_national_tract_rates.R since
# it's pure attribute lookup (no spatial joins) and much faster to rebuild.
#
# Sources, all either bundled with tigris or already cached locally by the
# national-rates build:
#   - county:  tigris::fips_codes (bundled, no network call)
#   - puma:    tigris::pumas() per state (tigris-cached after first fetch)
#   - cd:      the CD cartographic boundary file already downloaded to
#              data/cache/cd_boundaries/ by build_national_tract_rates.R
#
# Run once locally: Rscript code/build_geo_names.R

library(dplyr)
library(sf)
library(purrr)
library(stringr)

options(tigris_use_cache = TRUE)

root_dir  <- "."
out_path  <- file.path(root_dir, "data/cache/geo_names.rds")

# ---- county ----
data("fips_codes", package = "tigris")
county_names <- fips_codes |>
  dplyr::transmute(
    GEOID = paste0(state_code, county_code),
    name  = paste0(county, ", ", state)
  )

# ---- congressional district ----
cd_zip_path <- file.path(root_dir, "data/cache/cd_boundaries/cb_2022_us_cd118_500k.zip")
stopifnot(file.exists(cd_zip_path))

state_names <- fips_codes |> dplyr::distinct(state_code, state_name)

cd_names <- sf::st_read(paste0("/vsizip/", cd_zip_path), quiet = TRUE) |>
  sf::st_drop_geometry() |>
  dplyr::left_join(state_names, by = c("STATEFP" = "state_code")) |>
  dplyr::transmute(
    GEOID = GEOID,
    name  = paste0(state_name, " ", NAMELSAD)
  )

# ---- PUMA ----
# No bundled national PUMA name table; NAME10 on each state's own PUMA
# cartographic boundary fetch is the readable name (e.g. "Detroit City
# (Northeast)"). tigris caches each state's fetch locally after first use.
state_abbrs <- fips_codes |> dplyr::distinct(state) |> dplyr::pull(state)
state_abbrs <- state_abbrs[nchar(state_abbrs) == 2]  # drop territories without PUMAs

message("Fetching PUMA names for ", length(state_abbrs), " states...")

puma_names <- purrr::map_dfr(state_abbrs, function(st) {
  message("  ", st)
  tryCatch(
    tigris::pumas(state = st, year = 2019, cb = TRUE) |>
      sf::st_drop_geometry() |>
      dplyr::transmute(GEOID = GEOID10, name = NAME10),
    error = function(e) {
      message("    skipped (", conditionMessage(e), ")")
      tibble::tibble(GEOID = character(0), name = character(0))
    }
  )
})

geo_names <- list(county = county_names, puma = puma_names, cd = cd_names)

saveRDS(geo_names, out_path)
message("Saved: ", out_path, " (", round(file.size(out_path) / 1e6, 2), " MB)")
message("county: ", nrow(county_names), " | puma: ", nrow(puma_names), " | cd: ", nrow(cd_names))
