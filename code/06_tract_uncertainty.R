# code/06_tract_uncertainty.R
#
# B4, part 2. Tract-level uncertainty, done correctly.
#
# Two components, deliberately reported apart:
#
#   HPS. Push each of the 80 replicate weight sets through to the tract and
#     take the BRR variance THERE. Summing per-cell variances instead would
#     assume the ~38 age-sex cells have independent errors -- they don't, they
#     come from one sample -- and understates the tract CV by about sqrt(38).
#
#   ACS. The tract age-sex counts are themselves estimates. B01001 5-year MOEs
#     at this grain are large (median MOE/estimate ~0.7). Combined in quadrature
#     across cells, per Census guidance.
#
# Both remain a LOWER BOUND: the PUMA reweighting, IDW smoothing, capacity cap
# and Dirichlet composition are all unquantified here.

root_dir <- "."
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(stringr)
  library(purrr); library(tidycensus)
})
source("code/figures/theme_gaydar.R")
source("code/helpers.R")

R_REP <- 80; VAR_MULT <- 1/(R_REP * 0.25)
STATES <- c("CA","TX","NY","MS","WY","VT")
FIPS <- c(1,2,4,5,6,8,9,10,12,13,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,
          33,34,35,36,37,38,39,40,41,42,44,45,46,47,48,49,50,51,53,54,55,56)
abb2fips <- setNames(c(FIPS, 11L), c(state.abb, "DC"))

rep_wide <- as.data.table(readRDS("data/hps/rate_replicates.rds"))
repcols  <- paste0("PWEIGHT", 1:R_REP)

out <- map_dfr(STATES, function(st) {
  message("  ", st)
  acs <- get_acs(geography="tract", variables=B01001_VARS, state=st, year=2023,
                 survey="acs5", output="tidy", geometry=FALSE) |>
    filter(variable != "total") |>
    mutate(sex = case_when(str_starts(variable,"m_") ~ "M", str_starts(variable,"f_") ~ "F"),
           acs_bin = str_remove(variable, "^[mf]_")) |>
    filter(!is.na(sex), acs_bin %in% ADULT_BINS) |>
    transmute(GEOID, sex, acs_bin, pop = estimate, pop_se = moe/1.645)

  rr <- rep_wide[EST_ST == abb2fips[[st]]]
  keys <- rr[, .(sex, acs_bin)]
  Rmat <- as.matrix(rr[, ..repcols])          # cells x 80
  r0   <- rr$r0                                # cells

  # tract x cell population matrix, columns aligned to Rmat rows
  acs <- acs |> mutate(cell = paste(sex, acs_bin))
  cell_order <- paste(keys$sex, keys$acs_bin)
  P <- acs |> select(GEOID, cell, pop) |>
    pivot_wider(names_from = cell, values_from = pop, values_fill = 0)
  Pse <- acs |> select(GEOID, cell, pop_se) |>
    pivot_wider(names_from = cell, values_from = pop_se, values_fill = 0)
  keep <- cell_order[cell_order %in% names(P)]
  idx  <- match(keep, cell_order)
  Pm   <- as.matrix(P[, keep]); Psem <- as.matrix(Pse[, keep])

  l0  <- as.vector(Pm %*% r0[idx])                       # point estimate
  L   <- Pm %*% Rmat[idx, , drop = FALSE]                # tract x 80
  vhps<- VAR_MULT * rowSums((L - l0)^2)
  vacs<- as.vector(Psem^2 %*% (r0[idx]^2))               # quadrature over cells

  tibble(state_abbr = st, GEOID = P$GEOID, pop = rowSums(Pm),
         lgbt = l0, var_hps = vhps, var_acs = vacs)
}) |>
  filter(pop > 0, lgbt > 0) |>
  mutate(var_both = var_hps + var_acs,
         cv_hps  = sqrt(var_hps)/lgbt,
         cv_acs  = sqrt(var_acs)/lgbt,
         cv_both = sqrt(var_both)/lgbt,
         acs_share = var_acs/var_both)

cat(sprintf("\n=== B4: tract uncertainty, %s tracts (%s) ===\n",
            format(nrow(out), big.mark=","), paste(STATES, collapse="/")))
cat(sprintf("median CV -- HPS only %.3f | ACS only %.3f | combined %.3f  (%.0f%% wider than HPS alone)\n",
  median(out$cv_hps), median(out$cv_acs), median(out$cv_both),
  100*(median(out$cv_both)/median(out$cv_hps)-1)))
cat(sprintf("median share of variance from ACS: %.0f%%\n", 100*median(out$acs_share)))
cat("\nBy tract adult-population quintile:\n")
print(as.data.frame(out |> mutate(q = ntile(pop,5)) |> group_by(q) |>
  summarise(tracts=n(), med_pop=round(median(pop)),
            cv_hps=round(median(cv_hps),3), cv_acs=round(median(cv_acs),3),
            cv_both=round(median(cv_both),3),
            acs_pct=round(100*median(acs_share)), .groups="drop")))
saveRDS(out, "data/cache/b4_tract_uncertainty.rds")
cat("\nwrote data/cache/b4_tract_uncertainty.rds\n")
