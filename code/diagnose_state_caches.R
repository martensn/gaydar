# code/diagnose_state_caches.R
# One-off diagnostic: load every precomputed tract_state cache and check for
# red flags (invalid/empty geometry, NA blowups, coordinate anomalies that
# would indicate antimeridian issues, missing CBSA joins).
# Run from the gaydar project root: Rscript code/diagnose_state_caches.R

library(sf)
library(dplyr)

cache_dir <- "data/cache/tract_state"
files <- list.files(cache_dir, pattern = "^tract_.*\\.rds$", full.names = TRUE)
states <- gsub("tract_(.*)_2023.*", "\\1", basename(files))

results <- list()

for (i in seq_along(files)) {
  st <- states[i]
  f <- files[i]

  res <- tryCatch({
    obj <- readRDS(f)

    n <- nrow(obj)
    is_sf <- inherits(obj, "sf")
    geom_types <- if (is_sf) table(as.character(sf::st_geometry_type(obj))) else NULL
    n_invalid <- if (is_sf) sum(!sf::st_is_valid(obj), na.rm = TRUE) else NA
    n_empty <- if (is_sf) sum(sf::st_is_empty(obj)) else NA

    # numeric columns: check for NA / negative / non-finite blowups
    num_cols <- names(obj)[sapply(obj, is.numeric)]
    na_counts <- sapply(num_cols, function(cn) sum(is.na(obj[[cn]])))
    neg_counts <- sapply(num_cols, function(cn) sum(obj[[cn]] < 0, na.rm = TRUE))
    inf_counts <- sapply(num_cols, function(cn) sum(!is.finite(obj[[cn]]) & !is.na(obj[[cn]])))

    worst_na_col <- if (length(na_counts) > 0) names(which.max(na_counts)) else NA
    worst_na_n   <- if (length(na_counts) > 0) max(na_counts) else 0

    # CBSA coverage, if a cbsa-ish column exists
    cbsa_col <- grep("cbsa", names(obj), value = TRUE, ignore.case = TRUE)
    cbsa_na_pct <- if (length(cbsa_col) > 0) {
      round(100 * mean(is.na(obj[[cbsa_col[1]]])), 1)
    } else NA_real_

    # bounding box, to catch antimeridian-crossing weirdness (Alaska)
    bbox <- if (is_sf) sf::st_bbox(sf::st_transform(obj, 4326)) else NULL
    xrange <- if (!is.null(bbox)) bbox["xmax"] - bbox["xmin"] else NA

    list(
      state = st, ok = TRUE, n_tracts = n, is_sf = is_sf,
      geom_types = paste(names(geom_types), geom_types, sep = "=", collapse = ", "),
      n_invalid_geom = n_invalid, n_empty_geom = n_empty,
      worst_na_col = worst_na_col, worst_na_n = worst_na_n,
      any_negative = any(neg_counts > 0), any_nonfinite = any(inf_counts > 0),
      cbsa_col = if (length(cbsa_col) > 0) cbsa_col[1] else NA,
      cbsa_na_pct = cbsa_na_pct,
      lon_range = round(xrange, 1),
      error = NA_character_
    )
  }, error = function(e) {
    list(state = st, ok = FALSE, n_tracts = NA, is_sf = NA, geom_types = NA,
         n_invalid_geom = NA, n_empty_geom = NA, worst_na_col = NA, worst_na_n = NA,
         any_negative = NA, any_nonfinite = NA, cbsa_col = NA, cbsa_na_pct = NA,
         lon_range = NA, error = conditionMessage(e))
  })

  results[[st]] <- res
  cat(sprintf("[%2d/%2d] %s: %s\n", i, length(files), st,
              if (isTRUE(res$ok)) "loaded" else paste("ERROR:", res$error)))
}

df <- dplyr::bind_rows(results)
saveRDS(df, "data/cache/diagnostic_report.rds")
write.csv(df, "data/cache/diagnostic_report.csv", row.names = FALSE)

cat("\n=== FLAGGED STATES (any red flag) ===\n")
flagged <- df %>%
  filter(
    !ok |
    n_invalid_geom > 0 |
    n_empty_geom > 0 |
    worst_na_n > 0 |
    any_negative |
    any_nonfinite |
    (!is.na(cbsa_na_pct) & cbsa_na_pct > 50) |
    (!is.na(lon_range) & lon_range > 60)
  )
print(as.data.frame(flagged), row.names = FALSE)

cat("\n=== ALL 51 SUMMARY ===\n")
print(as.data.frame(df[, c("state","ok","n_tracts","n_invalid_geom","n_empty_geom",
                            "worst_na_col","worst_na_n","cbsa_na_pct","lon_range")]),
      row.names = FALSE)
