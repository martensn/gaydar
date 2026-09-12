# code/05_uncertainty_b4.R
#
# B4 -- uncertainty. The estimates currently carry none.
#
# TWO variance components, and it matters that they are reported separately:
#
#   (1) HPS sampling error in the state x sex x age identification rates.
#       Recoverable from the 80 replicate person weights (BRR, Fay k = 0.5, so
#       Var = 1/(R(1-k)^2) * sum (theta_r - theta_0)^2 = (1/20) * sum(...)).
#
#   (2) ACS sampling error in the tract age-sex counts the rates are projected
#       onto. B01001 5-year tract estimates carry their own MOEs, often large
#       for thin age-sex cells. Component (1) alone is a LOWER BOUND on tract
#       uncertainty -- quoting it as if it were the whole thing would understate
#       the interval, badly, in exactly the small tracts where it matters most.
#
# Not captured here, and worth stating in the paper: the PUMA reweighting, the
# IDW smoothing, the capacity cap and the Dirichlet composition all add further
# uncertainty. So even the combined figure remains a lower bound.

root_dir <- "."
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(stringr)
  library(purrr); library(tidycensus)
})
source("code/figures/theme_gaydar.R")
source("code/helpers.R")

R_REP <- 80; FAY_K <- 0.5; VAR_MULT <- 1 / (R_REP * (1 - FAY_K)^2)   # = 1/20

# ---------- (1) HPS replicate variance in state x sex x bin rates ----------
puf_files <- list.files("data/hps", pattern="^pulse[0-9]{4}_puf_[0-9]+\\.csv$", full.names=TRUE)
rep_files <- list.files("data/hps", pattern="^pulse[0-9]{4}_repwgt_puf_[0-9]+\\.csv$", full.names=TRUE)
repcols <- paste0("PWEIGHT", 1:R_REP)

message("loading ", length(puf_files), " PUF + replicate wave pairs ...")
cells <- rbindlist(map2(sort(puf_files), sort(rep_files), function(pf, rf) {
  p <- fread(pf, select=c("SCRAM","WEEK","EST_ST","EGENID_BIRTH","AGENID_BIRTH",
                          "TBIRTH_YEAR","SEXUAL_ORIENTATION","PWEIGHT"))
  r <- fread(rf, select=c("SCRAM","WEEK", repcols))
  d <- merge(p, r, by=c("SCRAM","WEEK"))
  d <- d[EGENID_BIRTH %in% 1:2 & PWEIGHT > 0]
  d[, sex := fifelse(EGENID_BIRTH==1,"M","F")]
  d[, age := 2023 - TBIRTH_YEAR]
  d[, y := as.integer(SEXUAL_ORIENTATION %in% c(1,3,4))]
  d[age >= 18]
}), fill=TRUE)

bins <- as.data.table(acs_age_bins)[, .(age_lo, age_hi, acs_bin)]
cells[, acs_bin := bins$acs_bin[findInterval(age, bins$age_lo)]]
cells <- cells[!is.na(acs_bin)]
message("pooled respondents: ", format(nrow(cells), big.mark=","))

pt  <- cells[, .(r0 = sum(PWEIGHT*y)/sum(PWEIGHT), n = .N), by=.(EST_ST, sex, acs_bin)]
rep <- cells[, lapply(.SD, function(w) sum(w*y)/sum(w)),
             by=.(EST_ST, sex, acs_bin), .SDcols=repcols]
rl  <- melt(rep, id.vars=c("EST_ST","sex","acs_bin"), value.name="r_rep")
v   <- merge(rl, pt, by=c("EST_ST","sex","acs_bin"))
rate_var <- v[, .(rate_se = sqrt(VAR_MULT * sum((r_rep - r0)^2))), by=.(EST_ST, sex, acs_bin)]
rate_var <- merge(rate_var, pt, by=c("EST_ST","sex","acs_bin"))
rate_var[, rate_cv := rate_se / pmax(r0, 1e-9)]

cat("\n=== (1) HPS replicate SEs on state x sex x age identification rate ===\n")
cat(sprintf("cells: %d   median rate %.3f%%   median SE %.3f pp   median CV %.3f\n",
    nrow(rate_var), 100*median(rate_var$r0), 100*median(rate_var$rate_se),
    median(rate_var$rate_cv)))
print(round(quantile(rate_var$rate_cv, c(.1,.25,.5,.75,.9,.99)), 3))
cat("\nworst-identified cells (highest CV):\n")
print(as.data.frame(rate_var[order(-rate_cv)][1:5,
  .(EST_ST, sex, acs_bin, n, rate = round(100*r0,2), se_pp = round(100*rate_se,2),
    cv = round(rate_cv,2))]))
saveRDS(rate_var, "data/hps/rate_replicate_se.rds")
cat("\nwrote data/hps/rate_replicate_se.rds\n")

# Save the per-replicate rate table too. Tract variance must be computed by
# pushing each replicate through to the tract and taking the BRR variance
# there -- NOT by summing per-cell variances. All cells in a state come from
# one sample, so their errors are strongly positively correlated; treating
# them as independent lets errors cancel and shrinks the tract CV by roughly
# sqrt(n_cells) (~6x here). Replicate propagation carries the covariance for
# free. See code/06_tract_uncertainty.R.
rep_wide <- merge(rep, pt[, .(EST_ST, sex, acs_bin, r0)], by=c("EST_ST","sex","acs_bin"))
saveRDS(rep_wide, "data/hps/rate_replicates.rds")
cat("wrote data/hps/rate_replicates.rds (", nrow(rep_wide), " cells x ", R_REP, " replicates )\n", sep="")
