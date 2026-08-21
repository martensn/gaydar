# code/build_national_geo_layers.R
#
# One-time offline build: dissolves the already-cached per-state tract
# layers up to county/PUMA/Congressional-District nationally, keeping
# geometry, so the map can render a single national choropleth at these
# geography levels instead of only the geocoded state.
#
# Tract level intentionally stays per-state at runtime -- rendering all
# ~85k tracts nationally as interactive Leaflet polygons would be far too
# slow. County (3,144), PUMA (2,351), and CD (436) are each comparable to
# or smaller than a single state's tract count the app already renders
# fine today, so those three go national.
#
# No population re-estimation happens here -- this is a pure geometric
# dissolve of numbers already computed, using the same group-by-and-sum
# logic as app.R's dissolve_geo(). Run once locally:
#   Rscript code/build_national_geo_layers.R

library(dplyr)
library(sf)
library(stringr)

root_dir  <- "."
cache_dir <- file.path(root_dir, "data/cache/tract_state")
out_path  <- file.path(root_dir, "data/cache/national_geo_layers.rds")

files <- list.files(cache_dir, pattern = "^tract_.*\\.rds$", full.names = TRUE)
stopifnot(length(files) > 0)

national_tract_rates <- readRDS(file.path(root_dir, "data/cache/national_tract_rates.rds"))
cd_lookup <- national_tract_rates |> dplyr::distinct(GEOID, cd_id)

# Same dissolve as app.R's dissolve_geo(): sum numeric columns within each
# group (sf's dplyr method auto-dissolves geometry), keep cbsa_id via
# first() since metro_bg() needs it and it's categorical not summable.
dissolve_geo <- function(sf_obj, group_col) {
  sf_obj |>
    dplyr::group_by(.data[[group_col]]) |>
    dplyr::summarise(
      dplyr::across(dplyr::where(is.numeric), \(x) sum(x, na.rm = TRUE)),
      cbsa_id = dplyr::first(cbsa_id),
      .groups = "drop"
    ) |>
    dplyr::rename(GEOID = dplyr::all_of(group_col)) |>
    dplyr::filter(!is.na(GEOID))
}

message("Building national county/PUMA/CD layers from ", length(files), " state caches...")

county_parts <- list()
puma_parts   <- list()
cd_parts     <- list()

for (f in files) {
  state_abbr <- stringr::str_match(basename(f), "^tract_([A-Z]{2})_")[, 2]
  message("  ", state_abbr)

  tract_sf <- readRDS(f)
  tract_sf <- tract_sf[sf::st_geometry_type(tract_sf) %in% c("POLYGON", "MULTIPOLYGON"), ]
  tract_sf <- sf::st_make_valid(tract_sf)
  tract_sf <- tract_sf[sf::st_is_valid(tract_sf) & !sf::st_is_empty(tract_sf), ]

  county_parts[[state_abbr]] <- dissolve_geo(tract_sf, "county_fips")
  puma_parts[[state_abbr]]   <- dissolve_geo(tract_sf, "puma_id")

  tract_cd <- tract_sf |> dplyr::left_join(cd_lookup, by = "GEOID")
  cd_parts[[state_abbr]] <- dissolve_geo(tract_cd, "cd_id")
}

message("Combining national layers...")
# county_fips/puma_id/cd_id are all state-FIPS-prefixed, so IDs are unique
# across states already -- a plain row-bind is safe, no re-grouping needed.
national_county <- dplyr::bind_rows(county_parts)
national_puma   <- dplyr::bind_rows(puma_parts)
national_cd     <- dplyr::bind_rows(cd_parts)

for (nm in c("national_county", "national_puma", "national_cd")) {
  obj <- get(nm)
  message(sprintf("  %-16s %5d rows, %s", nm, nrow(obj), format(object.size(obj), units = "MB")))
}

national_geo_layers <- list(county = national_county, puma = national_puma, cd = national_cd)
saveRDS(national_geo_layers, out_path)
message("Saved: ", out_path, " (", round(file.size(out_path) / 1e6, 2), " MB)")
