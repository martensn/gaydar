# code/build_all_state_caches.R
# One-time local batch job: precompute the tract-level expected-LGBTQ-population
# layer for every state + DC, so the deployed Shiny app never has to run the
# live ACS/tigris/rmapshaper pipeline on a user's request (was taking 8+ minutes
# per state on shinyapps.io and risked timeouts/OOM on larger states).
#
# Run from the gaydar project root: Rscript code/build_all_state_caches.R

root_dir <- "."

library(sf)
library(dplyr)
library(tidyr)
library(stringr)
library(tidycensus)
library(readr)
library(rmapshaper)
library(parallel)

options(tigris_use_cache = TRUE)

source("code/helpers.R")

rates <- readRDS(file.path(root_dir, "data/hps/hps_acs_rates.rds"))

cache_dir <- file.path(root_dir, "data/cache", "tract_state")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(root_dir, "data/cache", "build_log.txt")
write(paste0("=== Batch build started ", Sys.time(), " ==="), log_file, append = FALSE)

states <- c(state.abb, "DC")

build_one <- function(state) {
  cache_key  <- paste0("tract_", state, "_2023_calTRUE_g0.5.rds")
  cache_path <- file.path(cache_dir, cache_key)

  if (file.exists(cache_path)) {
    write(paste0(Sys.time(), " SKIP (already built): ", state), log_file, append = TRUE)
    return(invisible(NULL))
  }

  result <- tryCatch({
    t0 <- Sys.time()

    tmp <- build_tract_expected_layer(
      state_abbr      = state,
      year            = 2023,
      rates           = rates,
      use_calibration = TRUE,
      gamma           = 0.5
    )
    tmp <- tmp[sf::st_geometry_type(tmp) %in% c("POLYGON", "MULTIPOLYGON"), ]
    tmp <- rmapshaper::ms_simplify(tmp, keep = 0.05, keep_shapes = TRUE)
    tmp <- sf::st_make_valid(tmp)
    tmp <- tmp[!sf::st_is_empty(tmp), ]

    # match the extra validation app.R's bg_sf() reactive applies after loading
    tmp <- sf::st_make_valid(tmp)
    tmp <- tmp[sf::st_is_valid(tmp) & !sf::st_is_empty(tmp), ]

    saveRDS(tmp, cache_path)

    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    write(paste0(Sys.time(), " OK (", elapsed, "s, ", nrow(tmp), " tracts): ", state), log_file, append = TRUE)
    "ok"
  }, error = function(e) {
    write(paste0(Sys.time(), " ERROR: ", state, " -- ", conditionMessage(e)), log_file, append = TRUE)
    "error"
  })

  invisible(result)
}

# 4 parallel workers: balances wall-clock time against local memory pressure
# and Census API load. tigris_use_cache writes to distinct per-geography files
# so concurrent workers fetching different states' shapefiles don't collide.
invisible(parallel::mclapply(states, build_one, mc.cores = 4, mc.preschedule = FALSE))

write(paste0("=== Batch build finished ", Sys.time(), " ==="), log_file, append = TRUE)
