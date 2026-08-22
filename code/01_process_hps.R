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


# Empirical Bayes smoothing with age-local priors
# User can control the size of neighborhood
fit_dirichlet_mom <- function(P, w = NULL, eps = 1e-8) {
  
  P <- pmax(P, eps)
  P <- P / rowSums(P)
  
  if (is.null(w)) w <- rep(1, nrow(P))
  w_norm <- w / sum(w)
  
  m <- colSums(w_norm * P)
  
  v <- colSums(w_norm * (P - matrix(m, nrow(P), length(m), byrow = TRUE))^2)
  
  alpha0 <- mean(m * (1 - m) / pmax(v, eps) - 1)
  alpha0 <- max(alpha0, 1e-3)
  
  alpha <- m * alpha0
  
  list(alpha = alpha, alpha0 = alpha0)
}


smooth_age_dirichlet_decay <- function(df, tau = 3) {
  
  df %>%
    group_by(lgbt_cat, sex) %>%
    group_modify(function(d_group, keys) {
      
      ages <- sort(unique(d_group$age))
      
      purrr::map_dfr(ages, function(a) {
        
        d_focal <- d_group %>% filter(age == a)
        
        d_neighbors <- d_group %>%
          filter(age != a, n_eff > 0)
        
        if (nrow(d_neighbors) == 0) {
          return(d_focal %>% mutate(p_cat_eb = p_cat))
        }
        
        # apply smooth decay weights
        d_neighbors <- d_neighbors %>%
          mutate(
            w_adj = n_eff * exp(-abs(age - a) / tau)
          )
        
        # reshape to wide for Dirichlet estimation
        P_neighbors <- d_neighbors %>%
          select(age, gen_cat, p_cat, w_adj) %>%
          tidyr::pivot_wider(
            names_from  = gen_cat,
            values_from = p_cat
          )
        
        w <- P_neighbors$w_adj
        P_mat <- as.matrix(P_neighbors %>% select(-age, -w_adj))
        
        # estimate Dirichlet prior
        dir_fit <- fit_dirichlet_mom(P_mat, w)
        
        alpha <- dir_fit$alpha
        alpha0 <- sum(alpha)
        
        # shrink focal
        d_focal %>%
          arrange(gen_cat) %>%
          mutate(
            p_cat_eb =
              (n_eff * p_cat + alpha[match(gen_cat, names(alpha))]) /
              (n_eff + alpha0)
          )
      })
    }) %>%
    ungroup() %>%
    select(lgbt_cat,sex,age,gen_cat,p_cat_eb,n_eff,w_sum,w2_sum) %>%
    pivot_wider(names_from = "gen_cat", values_from = "p_cat_eb")
}


# Empirical Bayes smoothing of a scalar proportion (e.g. married share) with
# age-local priors, same decay-weighted-neighbor logic as
# smooth_age_dirichlet_decay() above but for a single Beta share instead of a
# Dirichlet-distributed composition vector.
smooth_age_beta_decay <- function(df, tau = 3) {

  df %>%
    group_by(lgbt_cat, sex, gen_cat) %>%
    group_modify(function(d_group, keys) {

      ages <- sort(unique(d_group$age))

      purrr::map_dfr(ages, function(a) {

        d_focal <- d_group %>% filter(age == a)

        d_neighbors <- d_group %>%
          filter(age != a, n_eff > 0)

        if (nrow(d_neighbors) == 0) {
          return(d_focal %>% mutate(married_shr_eb = married_shr))
        }

        d_neighbors <- d_neighbors %>%
          mutate(w_adj = n_eff * exp(-abs(age - a) / tau))

        fb <- fit_beta_mom(d_neighbors$married_shr, w = d_neighbors$w_adj)

        d_focal %>%
          mutate(
            married_shr_eb = (n_eff * married_shr + fb$alpha) / (n_eff + fb$alpha + fb$beta)
          )
      })
    }) %>%
    ungroup()
}


file_list <- list.files(file.path(root_dir,"data/hps"))
file_list <- file_list[str_detect(file_list, ".csv") & !str_detect(file_list, "repwgt")]

acs_age_bins <- tibble::tribble(
  ~age_lo, ~age_hi, ~acs_bin, ~acs_bin_col,
  18,  19, "18_19", "18_19",
  20,  20, "20", "20",
  21,  21, "21", "21",
  22,  24, "22_24", "22_24",
  25,  29, "25_29", "25_29",
  30,  34, "30_34", "30_34",
  35,  39, "35_39", "35_39",
  40,  44, "40_44", "40_44",
  45,  49, "45_49", "45_49",
  50,  54, "50_54", "50_54",
  55,  59, "55_59", "55_59",
  60,  61, "60_61", "60_61",
  62,  64, "62_64", "62+",
  65,  66, "65_66", "62+",
  67,  69, "67_69", "62+",
  70,  74, "70_74", "62+",
  75,  79, "75_79", "62+",
  80,  84, "80_84", "62+",
  85, 120, "85+", "62+"
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
    gender_cat = case_when(
      GENID_DESCRIBE == 3 & EGENID_BIRTH == 1 ~ "Trans woman",
      GENID_DESCRIBE == 3 & EGENID_BIRTH == 2 ~ "Trans man",
      GENID_DESCRIBE == 2 & EGENID_BIRTH == 1 ~ "Trans woman",
      GENID_DESCRIBE == 1 & EGENID_BIRTH == 2 ~ "Trans man",
      GENID_DESCRIBE == 1 & EGENID_BIRTH == 1 ~ "Cis man",
      GENID_DESCRIBE == 2 & EGENID_BIRTH == 2 ~ "Cis woman",
      GENID_DESCRIBE == 4 & EGENID_BIRTH == 1 ~ "NB AMAB",
      GENID_DESCRIBE == 4 & EGENID_BIRTH == 2 ~ "NB AFAB"
    ),
    # SEXUAL_ORIENTATION: 1 Gay/lesbian, 2 Straight, 3 Bisexual,
    # 4 Something else, 5 I don't know. Category 5 ("I don't know") is
    # treated as not-confidently-LGBT rather than folded into "Queer" --
    # someone unsure enough to answer "I don't know" isn't a confident LGBT
    # identification, and merging it into Queer produced an implausible
    # rising-with-age pattern (older respondents answering "I don't know"
    # for non-identity reasons swamping genuine "something else" identity).
    lgbt_cat = case_when(
      #gender_cat %in% c("Trans man", "Trans woman") ~ "Trans",
      SEXUAL_ORIENTATION == 3 ~ "Bisexual",
      SEXUAL_ORIENTATION == 1 ~ "LG",
      SEXUAL_ORIENTATION == 4 ~ "Queer",
      TRUE ~ "Straight"
    ),
    y_lgbt = as.integer(lgbt_cat != "Straight"),
    # MS: 1 Married, 2 Widowed, 3 Divorced, 4 Separated, 5 Never married, -99 missing.
    # HPS has no unmarried-cohabiting-partner flag, so "married" is the closest
    # available proxy for "coupled" -- used as the PUMA-reweighting blend weight.
    married = case_when(
      MS == 1          ~ 1L,
      MS %in% 2:5       ~ 0L,
      TRUE              ~ NA_integer_
    )
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

gen_sex_flow = hps_micro %>% 
  # Top code age at 62
  # I'd like to top code at a higher age, but the Dirichlet smoothing won't behave
  # because we start encountering missing data for the 62 year old cohort
  mutate(age = if_else(age > 62,62,age)) %>%
  #left_join(collapsed_age_bins, by = join_by(age >= age_lo, age <= age_hi)) %>%
  group_by(lgbt_cat,gender_cat,sex,age) %>% 
  summarize(n = sum(PWEIGHT)) %>%  
  filter(!is.na(gender_cat)) %>% 
  group_by(lgbt_cat,sex,age) %>% 
  mutate(p_cat = n/sum(n),
         gen_cat = case_when(gender_cat %in% c("Cis man","Cis woman") ~ "cis",
                             gender_cat %in% c("NB AFAB", "NB AMAB") ~ "non_binary",
                             gender_cat %in% c("Trans man", "Trans woman") ~ "trans"),
         w_sum  = sum(n, na.rm = TRUE),
         w2_sum = sum(n^2, na.rm = TRUE),
         n_eff  = ifelse(w2_sum > 0, (w_sum^2) / w2_sum, NA_real_))
         

# Tau is increasing in smoothness
gsf_eb <- smooth_age_dirichlet_decay(gen_sex_flow, tau = 3)
saveRDS(gsf_eb,file.path(root_dir,"data/hps","age_gender_sexuality.rds"))

# Merge on ACS age bins and take weighte averages to compute binned shares
gsf_eb_bin = gsf_eb %>%
  left_join(acs_age_bins, by = join_by(age >= age_lo, age <= age_hi)) %>%
  group_by(lgbt_cat,sex,acs_bin_col) %>%
  summarize(cis = weighted.mean(cis, w=w_sum),
            non_binary = weighted.mean(non_binary, w=w_sum),
            trans = weighted.mean(trans, w=w_sum)) %>%
  select(lgbt_cat,sex,acs_bin_col,cis,non_binary,trans)

saveRDS(gsf_eb_bin,file.path(root_dir,"data/hps","age_bin_gender_sexuality.rds"))


# --- Nationwide married-share by (sex, gender identity, sexuality, age bin) ---
# Used by the PUMA-reweighting step (code/helpers.R) as the per-cell blend
# weight between the same-sex-couples spatial signal and the singles spatial
# signal. Nationwide (not state-varying) -- HPS cells are already thin once
# split by sex/gender/sexuality/age, splitting further by state would leave
# too little sample to be usable.
married_flow <- hps_micro %>%
  # Same age topcode as gen_sex_flow above, for the same reason (Dirichlet/
  # Beta smoothing misbehaves once per-age sample gets too thin near 62+).
  mutate(age = if_else(age > 62, 62, age)) %>%
  filter(!is.na(gender_cat), !is.na(married)) %>%
  mutate(
    gen_cat = case_when(
      gender_cat %in% c("Cis man","Cis woman") ~ "cis",
      gender_cat %in% c("NB AFAB", "NB AMAB")   ~ "non_binary",
      gender_cat %in% c("Trans man","Trans woman") ~ "trans"
    )
  ) %>%
  group_by(lgbt_cat, sex, gen_cat, age) %>%
  summarize(
    married_shr = weighted.mean(married, PWEIGHT, na.rm = TRUE),
    w_sum  = sum(PWEIGHT, na.rm = TRUE),
    w2_sum = sum(PWEIGHT^2, na.rm = TRUE),
    n_eff  = ifelse(w2_sum > 0, (w_sum^2) / w2_sum, NA_real_),
    .groups = "drop"
  )

# Tau matches gen_sex_flow's smoothing above
married_eb <- smooth_age_beta_decay(married_flow, tau = 3)

# Collapse to ACS age bins (62+ collapsed), same grain as age_bin_gender_sexuality.rds
married_share_by_cell <- married_eb %>%
  left_join(acs_age_bins, by = join_by(age >= age_lo, age <= age_hi)) %>%
  group_by(lgbt_cat, sex, gen_cat, acs_bin_col) %>%
  summarize(
    married_shr = weighted.mean(married_shr_eb, w = w_sum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Floor/ceiling so no cell relies entirely on one spatial signal
  mutate(married_shr = pmin(pmax(married_shr, 0.05), 0.95))

saveRDS(married_share_by_cell, file.path(root_dir,"data/hps","married_share_by_cell.rds"))


# Compute relative shares of LGBT sub-populations by age-sex, pooled
# nationally rather than by state. The overall LGBT identification RATE
# (lgbt_shr, below) is the primary geographic signal and stays state-level,
# Gallup-calibrated -- but the *mix* of LG vs Bisexual vs Queer among people
# who do identify as LGBT is a much thinner slice of the HPS sample once
# split by sex x state x age bin x category, and cross-state Beta-shrinkage
# wasn't enough to make individual states' composition estimates trustworthy.
# Pooling nationally trades away state-level composition variation (which
# was mostly noise) for a much larger effective sample per age bin.
hps_micro_comp <- hps_micro %>%
  left_join(acs_age_bins, by = join_by(age >= age_lo, age <= age_hi)) %>%
  filter(y_lgbt == 1)

comp_week <- hps_micro_comp %>%
  group_by(sex, acs_bin, WEEK, lgbt_cat) %>%
  summarise(
    w_cat = sum(PWEIGHT, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(sex, acs_bin, WEEK) %>%
  mutate(
    w_total = sum(w_cat),
    shr_cat = w_cat / pmax(w_total, 1e-12)
  )

# Conditional probabilities of LGBT sub-populations within a given sex-age cell
comp_rates <- comp_week %>%
  group_by(sex, acs_bin, lgbt_cat) %>%
  summarise(
    shr_cat = mean(shr_cat, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  transmute(sex, acs_bin, lgbt_cat, shr_cat)

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

# No cross-state shrinkage needed any more -- comp_rates is already national.
# Just renormalize so LG/Bisexual/Queer/Straight sum to 1 within each
# (sex, acs_bin) cell (week-averaging above isn't guaranteed to land exactly
# on 1, since not every category necessarily appears in every week).
comp_rates_eb <- comp_rates %>%
  filter(!is.na(lgbt_cat), !is.na(shr_cat)) %>%
  group_by(sex, acs_bin) %>%
  mutate(
    shr_cat_eb = shr_cat / sum(shr_cat, na.rm = TRUE)
  ) %>%
  ungroup()

# Unconditional rates of identity by population. lgbt_shr (overall LGBT
# identification rate) stays state-level; shr_cat_eb (the LG/Bisexual/Queer
# mix conditional on identifying as LGBT) is now the same nationwide
# composition applied to every state's own identification rate.
rates_final <- hps_gallup_rates %>%
  left_join(
    comp_rates_eb %>% select(sex, acs_bin, lgbt_cat, shr_cat_eb),
    by = c("sex","acs_bin")
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










