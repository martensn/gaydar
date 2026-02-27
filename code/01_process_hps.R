library(data.table)
library(sf)
library(tidycensus)
library(stringr)
library(dplyr)
library(tidyr)

root_dir = "/Users/nickmartens/Desktop/gaydar"
#root_dir = "/Users/nickmartens/Library/Mobile Documents/com~apple~CloudDocs/Desktop/gaydar"
ref_year <- 2023

# Function used to compute parameters for Empirical Bayes shrinkage
fit_beta_mom <- function(p, w = NULL, eps = 1e-6) {
  p <- pmin(pmax(p, eps), 1 - eps)
  
  if (is.null(w)) {
    w <- rep(1, length(p))
  }
  
  # weighted mean
  m <- weighted.mean(p, w, na.rm = TRUE)
  
  # weighted variance
  w_norm <- w / sum(w, na.rm = TRUE)
  v <- sum(w_norm * (p - m)^2, na.rm = TRUE)
  
  # guard against degeneracy
  v <- max(v, 1e-8)
  v_max <- m * (1 - m)
  if (v >= v_max) v <- 0.99 * v_max
  
  k <- m * (1 - m) / v - 1
  alpha <- m * k
  beta  <- (1 - m) * k
  
  list(alpha = alpha, beta = beta, mean = m, var = v)
}

gen_from_age_2023 <- function(age) {
  dplyr::case_when(
    age <= 26 ~ "GenZ",
    age <= 42 ~ "Millennial",
    age <= 58 ~ "GenX",
    age <= 77 ~ "Boomer",
    TRUE      ~ "Silent"
  )
}

invlogit <- function(x) 1 / (1 + exp(-x))
logit <- function(p) log(p / (1 - p))

soft_calibrate_probs <- function(df, p0, w, X, margins, lambda = 1e4, maxit = 400) {
  eta0 <- logit(pmin(pmax(p0, 1e-6), 1 - 1e-6))
  
  obj <- function(theta) {
    eta <- eta0 + as.numeric(X %*% theta)
    p   <- invlogit(eta)
    
    # keep close to baseline (weighted L2 on logit shifts)
    #move <- sum(w * (eta - eta0)^2)
    move <- sum(w * (eta - eta0)^2) / sum(w)
    
    # margin penalty
    pen <- 0
    for (m in margins) {
      idx <- m$idx
      t   <- m$target
      mw  <- if (!is.null(m$weight)) m$weight else 1
      phat <- sum(w[idx] * p[idx]) / sum(w[idx])
      pen <- pen + mw * (phat - t)^2
    }
    move + lambda * pen
  }
  
  theta0 <- rep(0, ncol(X))
  fit <- optim(theta0, obj, method = "BFGS", control = list(maxit = maxit, reltol = 1e-12))
  
  theta <- fit$par
  eta   <- eta0 + as.numeric(X %*% theta)
  p_adj <- invlogit(eta)
  
  list(theta = theta, p_adj = p_adj, fit = fit)
}

margin_errors <- function(p, w, margins) {
  sapply(margins, function(m) {
    idx <- m$idx
    phat <- sum(w[idx] * p[idx]) / sum(w[idx])
    phat - m$target
  })
}


file_list <- list.files(file.path(root_dir,"data/hps"))
file_list <- file_list[str_detect(file_list, ".csv") & !str_detect(file_list, "repwgt")]

acs_age_bins <- tibble::tribble(
  ~age_lo, ~age_hi, ~acs_bin,
  18,  19, "18_19",
  20,  20, "20",
  21,  21, "21",
  22,  24, "22_24",
  25,  29, "25_29",
  30,  34, "30_34",
  35,  39, "35_39",
  40,  44, "40_44",
  45,  49, "45_49",
  50,  54, "50_54",
  55,  59, "55_59",
  60,  61, "60_61",
  62,  64, "62_64",
  65,  66, "65_66",
  67,  69, "67_69",
  70,  74, "70_74",
  75,  79, "75_79",
  80,  84, "80_84",
  85, 120, "85+"
)

# Import state crosswalk from tigris
state_crosswalk <- tigris::states(cb = TRUE, year = 2023) |>
  st_drop_geometry() |>
  transmute(
    state_fips = as.numeric(STATEFP),
    state_abbr = STUSPS,
    state = NAME
  ) |>
  arrange(state_fips) 


# bind together
hps <- file.path(root_dir, "data/hps", file_list) |>
  lapply(data.table::fread) |>
  dplyr::bind_rows()


#SEXUAL_ORIENTATION == 3 ~ "Bisexual",
#SEXUAL_ORIENTATION == 4 ~ "Other Queer",
#SEXUAL_ORIENTATION == 5 ~ "Questioning"),
#Gender = case_when(GENID_DESCRIBE == 1 ~ "AMAB",
#GENID_DESCRIBE == 2 ~ "AFAB",
#GENID_DESCRIBE == 3 ~ "Trans",
#GENID_DESCRIBE == 4 ~ "GNC"))
hps_micro <- hps %>%
  as_tibble() %>%
  mutate(
    sex = case_when(
      EGENID_BIRTH == 1 ~ "M",
      EGENID_BIRTH == 2 ~ "F",
      TRUE ~ NA_character_
    ),
    age = ref_year - TBIRTH_YEAR,
    gen = gen_from_age_2023(age),
    
    lgbt_cat = case_when(
      GENID_DESCRIBE == 3 ~ "Trans",
      SEXUAL_ORIENTATION == 3     ~ "Bisexual",
      SEXUAL_ORIENTATION == 1     ~ "LG",
      SEXUAL_ORIENTATION %in% c(4,5) ~ "Queer",
      TRUE ~ "Straight"
    ),
    
    y_lgbt = as.integer(lgbt_cat != "Straight")
  ) %>%
  filter(!is.na(sex), !is.na(gen), !is.na(PWEIGHT), PWEIGHT > 0) %>%
  mutate(
    gen = factor(gen, levels = c("GenZ","Millennial","GenX","Boomer","Silent")),
    sex = factor(sex, levels = c("F","M")),
    lgbt_cat = factor(lgbt_cat, levels = c("LG","Bisexual","Queer","Trans","Straight"))
  )

# Soft targets (within 1 p.p.)
# https://news.gallup.com/poll/656708/lgbtq-identification-rises.aspx
targets <- list(
  overall = 0.093,
  gen = c(
    GenZ       = 0.231,
    Millennial = 0.141,
    GenX       = 0.051,
    Boomer     = 0.030,
    Silent     = 0.018
    # Boomer not explicitly in your snippets; add if you want by pulling from the chart
  ),
  gen_sex = list(
    GenZ       = c(F = 0.31, M = 0.12),
    Millennial = c(F = 0.18, M = 0.09)
  )
)

# Use gen main effects + gen:sex interactions (we’ll softly target only two gens’ sex splits)
X <- model.matrix(~ 0 + gen + gen:sex, data = hps_micro)
w <- hps_micro$PWEIGHT

# Base probabilities: start from raw indicator (smoothed)
p0 <- pmin(pmax(hps_micro$y_lgbt, 1e-4), 1 - 1e-4)


margins <- list(
  list(idx = rep(TRUE, nrow(hps_micro)), target = targets$overall, weight = 2),
  
  list(idx = hps_micro$gen == "GenZ",       target = targets$gen[["GenZ"]],       weight = 3),
  list(idx = hps_micro$gen == "Millennial", target = targets$gen[["Millennial"]], weight = 3),
  list(idx = hps_micro$gen == "GenX",       target = targets$gen[["GenX"]],       weight = 3),
  list(idx = hps_micro$gen == "Boomer",     target = targets$gen[["Boomer"]],       weight = 3),
  list(idx = hps_micro$gen == "Silent",     target = targets$gen[["Silent"]],     weight = 3),
  
  list(idx = hps_micro$gen == "GenZ" & hps_micro$sex == "F", target = targets$gen_sex$GenZ[["F"]], weight = 5),
  list(idx = hps_micro$gen == "GenZ" & hps_micro$sex == "M", target = targets$gen_sex$GenZ[["M"]], weight = 5),
  list(idx = hps_micro$gen == "Millennial" & hps_micro$sex == "F", target = targets$gen_sex$Millennial[["F"]], weight = 5),
  list(idx = hps_micro$gen == "Millennial" & hps_micro$sex == "M", target = targets$gen_sex$Millennial[["M"]], weight = 5)
)

for (lam in c(1e5)) {
  res <- soft_calibrate_probs(hps_micro, p0, w, X, margins, lambda = lam)
  err <- margin_errors(res$p_adj, w, margins)
  cat("lambda=", lam, " max|err|=", max(abs(err)), "\n")
}

best <- soft_calibrate_probs(hps_micro, p0, w, X, margins, lambda = 3e4)  # example
hps_micro$p_lgbt_adj <- best$p_adj

# Compute relative shares of LGBT sub-populations by age-sex-state
hps_micro_comp <- hps_micro %>%
  left_join(acs_age_bins, by = join_by(age >= age_lo, age <= age_hi)) %>%
  filter(y_lgbt == 1)

comp_week <- hps_micro_comp %>%
  group_by(sex, EST_ST, acs_bin, WEEK, lgbt_cat) %>%
  summarise(
    w_cat = sum(PWEIGHT, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(sex, EST_ST, acs_bin, WEEK) %>%
  mutate(
    w_total = sum(w_cat),
    shr_cat = w_cat / pmax(w_total, 1e-12)
  )

# Conditional probabilities of LGBT sub-populations with a given sex-state-age cell
comp_rates <- comp_week %>%
  left_join(state_crosswalk %>% select(state_fips,state_abbr), by = c("EST_ST"="state_fips")) %>%
  group_by(sex, state_abbr, acs_bin, lgbt_cat) %>%
  summarise(
    shr_cat = mean(shr_cat, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  transmute(sex, state_abbr, acs_bin, lgbt_cat, shr_cat)

# Compute sum of person weights and n_eff, which measures the quality of the signal
hps_gallup_week <- hps_micro %>%
  left_join(acs_age_bins, by = join_by(age >= age_lo, age <= age_hi)) %>%
  group_by(sex, EST_ST, acs_bin, WEEK) %>%
  summarise(
    # expected LGB “count” under calibrated probs
    lgbt_w = sum(PWEIGHT * p_lgbt_adj, na.rm = TRUE),
    denom_w = sum(PWEIGHT, na.rm = TRUE),
    w_sum  = sum(PWEIGHT, na.rm = TRUE),
    w2_sum = sum(PWEIGHT^2, na.rm = TRUE),
    n_eff  = ifelse(w2_sum > 0, (w_sum^2) / w2_sum, NA_real_),
    .groups = "drop"
  ) %>%
  mutate(lgbt_shr = lgbt_w / pmax(denom_w, 1e-12))

# Pool results from multiple weeks
hps_gallup_rates <- hps_gallup_week %>%
  group_by(sex, EST_ST, acs_bin) %>%
  summarise(
    n_wsum = mean(w_sum, na.rm = TRUE),
    n_eff  = mean(n_eff, na.rm = TRUE),
    lgbt_shr = mean(lgbt_shr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(state_crosswalk, by = c("EST_ST" = "state_fips")) %>%
  transmute(sex, state, state_abbr, acs_bin, n_wsum, n_eff, lgbt_shr)

# Empirical Bayes smoothing on the sub-population shares
# Shrinking towards averages by age across states 
comp_rates_eb <- comp_rates %>%
  left_join(
    hps_gallup_rates %>% select(sex, state_abbr, acs_bin, n_eff),
    by = c("sex","state_abbr","acs_bin")
  ) %>%
  filter(!is.na(lgbt_cat), !is.na(shr_cat), !is.na(n_eff), n_eff > 0) %>%
  group_by(sex, acs_bin, lgbt_cat) %>%
  group_modify(~{
    d <- .x
    # If only 1 state in this group, nothing to shrink
    if (nrow(d) < 2) return(d %>% mutate(shr_cat_eb = shr_cat))
    
    fb <- fit_beta_mom(d$shr_cat, w = d$n_eff)
    
    d %>% mutate(
      shr_cat_eb = (n_eff * shr_cat + fb$alpha) / (n_eff + fb$alpha + fb$beta)
    )
  }) %>%
  ungroup() %>%
  group_by(sex, state_abbr, acs_bin) %>%
  mutate(
    shr_cat_eb = shr_cat_eb / sum(shr_cat_eb, na.rm = TRUE)
  ) %>%
  ungroup()

# Unconditional rates of identity by population
rates_final <- hps_gallup_rates %>%
  left_join(
    comp_rates_eb %>% select(sex, state_abbr, acs_bin, lgbt_cat, shr_cat_eb),
    by = c("sex","state_abbr","acs_bin")
  ) %>%
  mutate(
    lgbt_cat_shr = lgbt_shr * shr_cat_eb
  ) %>%
  filter(!is.na(lgbt_cat)) %>%
  select(sex, state, state_abbr, acs_bin, lgbt_cat, lgbt_shr, lgbt_cat_shr, n_wsum) %>%
  pivot_wider(names_from = lgbt_cat,
              values_from = lgbt_cat_shr,
              values_fill = 0)

saveRDS(rates_final,file.path(root_dir,"data/hps","hps_acs_rates.rds"))
