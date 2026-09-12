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

# Optional CLI args, so an alternative build can be written alongside the live
# caches instead of over them:
#   Rscript code/build_all_state_caches.R [cache_subdir] [rates_rds]
# Both default to the current HPS behaviour, so a bare invocation is unchanged.
#   e.g. ... build_all_state_caches.R tract_state_agenid
#        ... build_all_state_caches.R tract_state_brfss data/brfss/brfss_acs_rates.rds
args <- commandArgs(trailingOnly = TRUE)
cache_subdir <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else "tract_state"
rates_file   <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else "data/hps/hps_acs_rates.rds"
# Optional 3rd arg: comma-separated state list, or a file containing one. The
# BRFSS arm covers only 41 states, and states with no rates would otherwise
# join to NA and write a cache full of NAs rather than failing loudly.
state_arg    <- if (length(args) >= 3 && nzchar(args[[3]])) args[[3]] else ""

rates <- readRDS(file.path(root_dir, rates_file))

cache_dir <- file.path(root_dir, "data/cache", cache_subdir)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  root_dir, "data/cache",
  if (cache_subdir == "tract_state") "build_log.txt"
  else paste0("build_log_", cache_subdir, ".txt")
)
cat("cache_dir:", cache_dir, "\nrates:", rates_file, "\nlog:", log_file, "\n")
write(paste0("=== Batch build started ", Sys.time(), " ==="), log_file, append = FALSE)

states <- c(state.abb, "DC")
if (nzchar(state_arg)) {
  states <- if (file.exists(state_arg)) strsplit(trimws(readLines(state_arg, warn = FALSE)[1]), ",")[[1]]
            else strsplit(state_arg, ",")[[1]]
  states <- trimws(states)
}
# Fail loudly rather than silently writing NA caches for uncovered states.
missing_rates <- setdiff(states, unique(rates$state_abbr))
if (length(missing_rates))
  stop("no rates for: ", paste(missing_rates, collapse = " "))
cat("building", length(states), "states\n")

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
