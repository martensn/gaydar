# R/helpers.R
# Utilities for building a tract-level sf layer with expected LGBTQ counts/shares.
library(sf)
library(dplyr)
library(tidyr)
library(stringr)
library(tidycensus)

CRS_LEAFLET <- 4326
root_dir <- "."
# ACS B01001 detailed age-sex counts
B01001_VARS <- c(
  total = "B01001_001",
  m_0_4  = "B01001_003", m_5_9  = "B01001_004", m_10_14="B01001_005",
  m_15_17="B01001_006", m_18_19="B01001_007", m_20  ="B01001_008",
  m_21  ="B01001_009", m_22_24="B01001_010", m_25_29="B01001_011",
  m_30_34="B01001_012", m_35_39="B01001_013", m_40_44="B01001_014",
  m_45_49="B01001_015", m_50_54="B01001_016", m_55_59="B01001_017",
  m_60_61="B01001_018", m_62_64="B01001_019", m_65_66="B01001_020",
  m_67_69="B01001_021", m_70_74="B01001_022", m_75_79="B01001_023",
  m_80_84="B01001_024", m_85up ="B01001_025",
  f_0_4  ="B01001_027", f_5_9  ="B01001_028", f_10_14="B01001_029",
  f_15_17="B01001_030", f_18_19="B01001_031", f_20  ="B01001_032",
  f_21  ="B01001_033", f_22_24="B01001_034", f_25_29="B01001_035",
  f_30_34="B01001_036", f_35_39="B01001_037", f_40_44="B01001_038",
  f_45_49="B01001_039", f_50_54="B01001_040", f_55_59="B01001_041",
  f_60_61="B01001_042", f_62_64="B01001_043", f_65_66="B01001_044",
  f_67_69="B01001_045", f_70_74="B01001_046", f_75_79="B01001_047",
  f_80_84="B01001_048", f_85up ="B01001_049"
)

# Adult ACS bins used in your pipeline
ADULT_BINS <- c(
  "18_19","20","21","22_24","25_29","30_34","35_39","40_44","45_49",
  "50_54","55_59","60_61","62_64","65_66","67_69","70_74","75_79","80_84","85+"
)

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

library(h3jsr)
library(rmapshaper)

# Aggregate a tract SF layer into H3 hexagons for choropleth display.
# res = 8 (~737m across) is appropriate for city-zoom; try res = 7 for a coarser view.
build_hex_layer <- function(tract_sf, res = 8) {
  tract_valid <- tract_sf |>
    dplyr::filter(!sf::st_is_empty(geometry)) |>
    sf::st_transform(4326) |>
    sf::st_make_valid() |>
    dplyr::filter(
      !sf::st_is_empty(geometry),
      sf::st_geometry_type(geometry) %in% c("POLYGON", "MULTIPOLYGON"),
      sf::st_is_valid(geometry)
    )

  centroids <- sf::st_centroid(tract_valid)
  h3_ids    <- h3jsr::point_to_cell(centroids, res = res)

  tract_valid |>
    sf::st_drop_geometry() |>
    dplyr::mutate(h3_id = h3_ids) |>
    dplyr::filter(!is.na(h3_id)) |>
    dplyr::group_by(h3_id) |>
    dplyr::summarise(dplyr::across(where(is.numeric), \(x) sum(x, na.rm = TRUE)),
                     .groups = "drop") |>
    dplyr::mutate(geometry = h3jsr::cell_to_polygon(h3_id, simple = FALSE)$geometry) |>
    sf::st_as_sf(crs = 4326)
}

get_state_choices <- function() {
  tigris::states(cb = TRUE, year = 2023) |>
    sf::st_drop_geometry() |>
    dplyr::arrange(NAME) |>
    dplyr::pull(STUSPS)
}



# Convert ACS wide B01001 into long with sex + acs_bin for adults only
acs_wide_to_long_adults <- function(tract_wide, tract_puma_xwalk) {
  long <- tract_wide |>
    sf::st_drop_geometry() |>
    inner_join(tract_puma_xwalk, by = "GEOID", relationship = "one-to-many") %>%
    tidyr::pivot_longer(cols = dplyr::all_of(names(B01001_VARS)),
                        names_to = "var", values_to = "pop") |>
    dplyr::mutate(
      pop = pop * w_puma,
      # parse sex + bin from var names like m_22_24, f_70_74
      sex = dplyr::case_when(
        stringr::str_starts(var, "m_") ~ "M",
        stringr::str_starts(var, "f_") ~ "F",
        TRUE ~ NA_character_
      ),
      acs_bin = stringr::str_remove(var, "^[mf]_"),
      # keep only adult bins for applying HPS rates
      is_adult_bin = acs_bin %in% ADULT_BINS
    ) |>
    dplyr::filter(is_adult_bin) |>
    dplyr::select(GEOID, puma_id, sex, acs_bin, pop, w_puma) |>
    dplyr::left_join(acs_age_bins %>% select(acs_bin, acs_bin_col), by = "acs_bin")

  long
}


get_state_couples_signal <- function(
    state_abbr,
    year,
    min_couples = 100,
    winsor = c(0.01, 0.99)
) {
  
  ss_vars <- c("B09019_011", "B09019_013")
  os_vars <- c("B09019_010", "B09019_012")
  
  acs <- tidycensus::get_acs(
    geography = "puma",
    variables = c(ss_vars, os_vars),
    state = state_abbr,
    year = 2019,
    survey = "acs5",
    geometry = FALSE,
    cache_table = TRUE
  )
  
  puma <- acs |>
    dplyr::mutate(
      puma_id = GEOID,
      ss = ifelse(variable %in% ss_vars, estimate, 0),
      os = ifelse(variable %in% os_vars, estimate, 0)
    ) |>
    dplyr::group_by(puma_id) |>
    dplyr::summarise(
      ss = sum(ss, na.rm = TRUE),
      os = sum(os, na.rm = TRUE),
      couples = ss + os,
      .groups = "drop"
    ) |>
    dplyr::filter(couples >= min_couples)
  
  if (nrow(puma) == 0) return(NULL)
  
  puma <- puma |>
    dplyr::mutate(p = ss / couples)
  
  qs <- stats::quantile(puma$p, probs = winsor, na.rm = TRUE)

  puma |>
    dplyr::mutate(
      p = pmin(pmax(p, qs[1]), qs[2]),
      puma_signal = qlogis(p),
      # Normalized PUMA-level share of same-sex couples within the state --
      # a probability distribution over PUMAs (sums to 1), used as the
      # "couples concentration" half of the cell-level blend in
      # build_tract_expected_layer().
      couple_share_norm = ss / sum(ss, na.rm = TRUE)
    ) |>
    dplyr::select(puma_id, puma_signal, couple_share_norm)
}


# ACS table B12002 ("Sex by Marital Status by Age") variable layout, verified
# against the 2022 ACS5 variables.json rather than assumed:
#
#   sex_start + 0            sex subtotal (unused)
#   sex_start + 1             "Never married" subtotal
#   sex_start + 2  .. +15     "Never married" x 14 age brackets
#   sex_start + 16            "Now married" subtotal
#   sex_start + 17            "Married, spouse present" subtotal
#   sex_start + 18 .. +31     "Married, spouse present" x 14 age brackets
#   sex_start + 32            "Married, spouse absent" subtotal
#   sex_start + 33            "Separated" subtotal
#   sex_start + 34 .. +47     "Separated" x 14 age brackets
#   sex_start + 48            "Other" (spouse absent, not separated) subtotal
#   sex_start + 49 .. +62     "Other" x 14 age brackets
#   sex_start + 63            "Widowed" subtotal
#   sex_start + 64 .. +77     "Widowed" x 14 age brackets
#   sex_start + 78            "Divorced" subtotal
#   sex_start + 79 .. +92     "Divorced" x 14 age brackets
#
# sex_start = 2 for Male, 95 for Female (each sex block is 93 variables).
# The 14 age brackets, in order: 15-17, 18-19, 20-24, 25-29, 30-34, 35-39,
# 40-44, 45-49, 50-54, 55-59, 60-64, 65-74, 75-84, 85+.
#
# "Married" here means "married, spouse present" specifically (i.e. actually
# cohabiting with a spouse) -- "married, spouse absent" is grouped with
# never-married/widowed/divorced into "not married" below, since those
# people aren't living with a partner either, which is the concept this
# signal is meant to approximate.
b12002_sex_start <- c(M = 2L, F = 95L)

b12002_age_brackets <- tibble::tibble(
  age_lo = c(15, 18, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 75, 85),
  age_hi = c(17, 19, 24, 29, 34, 39, 44, 49, 54, 59, 64, 74, 84, 120),
  offset = 0:13
)

# Crosswalk from B12002's native (coarser) age brackets to this pipeline's
# acs_bin_col. Where one B12002 bracket spans several acs_bin_col groups
# (20-24 covers "20"/"21"/"22_24"; 65-74/75-84/85+ all fold into "62+"),
# every matching acs_bin_col reuses the same B12002-derived PUMA-share curve
# -- B12002 can't distinguish them any finer, and this signal only needs to
# be a reasonable spatial proxy, not a precise headcount. The 60-64 bracket
# is assigned to "60_61" only (rather than split); "62+" is built from the
# 65-74/75-84/85+ brackets instead, mirroring how acs_age_bins already
# collapses those same ages into one "62+" bucket for the rest of the
# pipeline.
b12002_to_acs_bin_col <- tibble::tribble(
  ~age_lo, ~acs_bin_col,
  18, "18_19",
  20, "20",
  20, "21",
  20, "22_24",
  25, "25_29",
  30, "30_34",
  35, "35_39",
  40, "40_44",
  45, "45_49",
  50, "50_54",
  55, "55_59",
  60, "60_61",
  65, "62+",
  75, "62+",
  85, "62+"
)

get_state_singles_signal <- function(
    state_abbr,
    year,
    min_singles = 100
) {

  var_grid <- tidyr::expand_grid(
    sex = names(b12002_sex_start),
    b12002_age_brackets
  ) |>
    dplyr::mutate(
      sex_start   = b12002_sex_start[sex],
      never_var   = sprintf("B12002_%03d", sex_start + 2  + offset),
      married_var = sprintf("B12002_%03d", sex_start + 18 + offset),
      sep_var     = sprintf("B12002_%03d", sex_start + 34 + offset),
      other_var   = sprintf("B12002_%03d", sex_start + 49 + offset),
      widowed_var = sprintf("B12002_%03d", sex_start + 64 + offset),
      divorced_var= sprintf("B12002_%03d", sex_start + 79 + offset)
    ) |>
    dplyr::left_join(b12002_to_acs_bin_col, by = "age_lo", relationship = "many-to-many") |>
    dplyr::filter(!is.na(acs_bin_col))

  var_cols <- c("never_var","married_var","sep_var","other_var","widowed_var","divorced_var")
  all_vars <- unique(unlist(var_grid[var_cols]))

  acs <- tidycensus::get_acs(
    geography = "puma",
    variables = stats::setNames(all_vars, all_vars),
    state = state_abbr,
    year = 2019,
    survey = "acs5",
    geometry = FALSE,
    cache_table = TRUE
  )

  long <- var_grid |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(var_cols),
      names_to = "kind",
      values_to = "variable"
    ) |>
    dplyr::mutate(
      kind = dplyr::if_else(kind == "married_var", "married", "not_married_part")
    ) |>
    dplyr::select(sex, acs_bin_col, age_lo, kind, variable) |>
    dplyr::inner_join(acs |> dplyr::rename(puma_id = GEOID), by = "variable", relationship = "many-to-many")

  singles <- long |>
    # collapse to one row per (puma, sex, source B12002 age bracket, kind)
    # before folding brackets into acs_bin_col, so a bracket shared across
    # several acs_bin_col groups (e.g. 20-24) isn't summed multiple times
    dplyr::group_by(puma_id, sex, acs_bin_col, age_lo, kind) |>
    dplyr::summarise(estimate = sum(estimate, na.rm = TRUE), .groups = "drop") |>
    dplyr::group_by(puma_id, sex, acs_bin_col, kind) |>
    dplyr::summarise(estimate = sum(estimate, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = kind, values_from = estimate, values_fill = 0) |>
    dplyr::mutate(
      total = married + not_married_part,
      not_married = not_married_part
    ) |>
    dplyr::filter(not_married >= min_singles) |>
    dplyr::group_by(sex, acs_bin_col) |>
    dplyr::mutate(
      singles_share_norm = not_married / sum(not_married, na.rm = TRUE)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(puma_id, sex, acs_bin_col, singles_share_norm)

  if (nrow(singles) == 0) return(NULL)

  singles
}


# app.R uses this when selecting genders
metric_gender_cols <- list(
  lgbt = list(
    m  = "lgbt_m",
    w  = "lgbt_w",
    nb = "lgbt_nb"
  ),
  lg = list(
    m  = c("m_nt_lg", "m_t_lg"),
    w  = c("w_nt_lg", "w_t_lg"),
    nb = "nb_nt_lg"
  ),
  bi = list(
    m  = c("m_nt_bisexual", "m_t_bisexual"),
    w  = c("w_nt_bisexual", "w_t_bisexual"),
    nb = "nb_nt_bisexual"
  ),
  queer = list(
    m  = c("m_nt_queer", "m_t_queer"),
    w  = c("w_nt_queer", "w_t_queer"),
    nb = "nb_nt_queer"
  ),
  trans = list(
    m  = c("m_t_lg", "m_t_bisexual", "m_t_queer", "m_t_straight"),
    w  = c("w_t_lg", "w_t_bisexual", "w_t_queer", "w_t_straight"),
    nb = character(0)
  )
)



build_tract_expected_layer <- function(
    state_abbr,
    year,
    rates,
    use_calibration = TRUE,
    gamma = 1) {

  # --- Tract-level age-sex data ---
  acs_raw <- tidycensus::get_acs(
    geography = "tract",
    variables = B01001_VARS,
    state = state_abbr,
    year = year,
    survey = "acs5",
    geometry = TRUE,
    output = "tidy",
    cache_table = TRUE
  )
  
  tract_wide <- acs_raw |>
    dplyr::select(GEOID, variable, estimate, geometry) |>
    tidyr::pivot_wider(
      names_from = "variable",
      values_from = "estimate",
      values_fill = 0
    ) |>
    sf::st_as_sf()

  # --- Map tract → PUMA (spatial join once) ---
  pumas <- tigris::pumas(
    state = state_abbr,
    year = 2019,
    cb = TRUE
  ) |>
    dplyr::select(puma_id = GEOID10)
  
  cbsas_sf <- tigris::core_based_statistical_areas(
    year = 2019,
    cb = TRUE
  ) |>
    dplyr::select(cbsa_id = GEOID, cbsa_name = NAME) |>
    sf::st_make_valid()

  # A handful of states' CBSA polygons (seen for AZ, MN) contain a
  # degenerate ring (duplicate consecutive vertex) that st_make_valid()
  # doesn't fully repair under the s2 spherical backend, and later breaks
  # st_join(..., largest = TRUE)'s st_area() call below. GEOS's
  # st_make_valid() handles this class of degeneracy more robustly; use it
  # for repair, then switch back to s2 (relied on elsewhere in this
  # function for antimeridian-safe distance calculations).
  sf::sf_use_s2(FALSE)
  cbsas_sf <- sf::st_make_valid(cbsas_sf)
  sf::sf_use_s2(TRUE)

  cbsas <- cbsas_sf |> st_drop_geometry()

  county_cbsa_xwalk = read_csv(file.path(root_dir,"data/unified_cbsa.csv")) %>%
    select(cbsa_id = cbsa_code, county_fips = GeoFIPS) %>%
    left_join(cbsas, by = "cbsa_id")
  
  tract_clean <- tract_wide |>
    st_transform(st_crs(pumas)) |>
    sf::st_make_valid()

  pumas_clean <- pumas |>
    sf::st_make_valid()

  tract_puma_xwalk <- st_intersection(
    tract_clean |> select(GEOID),
    pumas_clean |> select(puma_id)
  ) |>
    # re-validate intersection
    sf::st_make_valid() |>              
    filter(!sf::st_is_empty(geometry)) |>
    mutate(area_int = as.numeric(sf::st_area(geometry))) |>
    group_by(GEOID) |>
    mutate(
      area_bg = sum(area_int),
      w_puma = area_int / area_bg
    ) |>
    ungroup() |>
    select(GEOID, puma_id, w_puma) |>
    mutate(county_fips = substr(GEOID,1,5)) |>
    left_join(county_cbsa_xwalk, by = "county_fips")

  # --- Spatial fallback for stale/missing county->CBSA crosswalk entries ---
  # unified_cbsa.csv is a static county-FIPS-keyed lookup and can go stale
  # when county-equivalent boundaries change (e.g. Connecticut replaced its
  # 8 counties with 9 planning regions in 2022; the crosswalk's rows for the
  # new region codes point at a placeholder "09999" CBSA that doesn't exist,
  # so cbsa_name comes back NA even though cbsa_id is technically non-NA).
  # Any tract without a real cbsa_name after the FIPS-based join is resolved
  # directly against the actual CBSA polygons, which is authoritative
  # regardless of whether the FIPS crosswalk has caught up for a given
  # county's current code.
  needs_fallback <- tract_puma_xwalk |>
    filter(is.na(cbsa_name)) |>
    distinct(GEOID) |>
    pull(GEOID)

  if (length(needs_fallback) > 0) {
    fallback_geoms <- tract_clean |>
      filter(GEOID %in% needs_fallback) |>
      select(GEOID) |>
      sf::st_transform(sf::st_crs(cbsas_sf))

    # largest = TRUE computes st_area() on the tract/CBSA intersection to
    # pick the biggest-overlap CBSA -- for some states (seen for AZ, MN)
    # that intersection produces a technically-degenerate ring (duplicate
    # vertex) even when both inputs are individually valid, which s2's
    # stricter area computation rejects outright. GEOS's area computation
    # tolerates it, so compute this one join under GEOS rather than s2
    # (relied on elsewhere in this function for antimeridian-safe distance
    # calculations, so only switched off for this call).
    sf::sf_use_s2(FALSE)
    fallback_match <- sf::st_join(
      fallback_geoms, cbsas_sf, join = sf::st_intersects, largest = TRUE
    ) |>
      sf::st_drop_geometry() |>
      distinct(GEOID, .keep_all = TRUE) |>
      rename(cbsa_id_fallback = cbsa_id, cbsa_name_fallback = cbsa_name)
    sf::sf_use_s2(TRUE)

    tract_puma_xwalk <- tract_puma_xwalk |>
      left_join(fallback_match, by = "GEOID") |>
      mutate(
        use_fallback = is.na(cbsa_name),
        cbsa_id   = dplyr::if_else(use_fallback, cbsa_id_fallback, cbsa_id),
        cbsa_name = dplyr::if_else(use_fallback, cbsa_name_fallback, cbsa_name)
      ) |>
      select(-cbsa_id_fallback, -cbsa_name_fallback, -use_fallback)
  }

  # --- Adult population long ---
  long <- acs_wide_to_long_adults(tract_wide, tract_puma_xwalk)
  
  rates <- readRDS(file.path(root_dir,"data/hps/hps_acs_rates.rds"))
  gsf_eb_bin <- readRDS(file.path(root_dir,"data/hps","age_bin_gender_sexuality.rds"))
  
  rates_state <- rates %>%
    filter(state_abbr == !!state_abbr) %>%
    select(-state_abbr)
  
  long2 <- long %>%
    # Remove age xsex x cbg cells with no population
    filter(pop > 0) %>%
    # Merge on state-age-sex estimates of sexuality
    left_join(rates_state, by = c("sex", "acs_bin")) %>%
    mutate(Straight = 1 - LG - Bisexual - Queer) %>%
    pivot_longer(cols = c(LG,Bisexual,Queer,Straight),
                 names_to = "lgbt_cat",
                 values_to = "lgbt_cat_shr") %>%
    left_join(gsf_eb_bin, by = c("acs_bin_col","lgbt_cat","sex")) %>%
    pivot_longer(cols = c(cis,non_binary,trans),
                 names_to = "gender_cat",
                 values_to = "gender_cat_shr") %>%
    mutate(gender = case_when(sex == "M" & gender_cat == "cis" ~ "m",
                              sex == "F" & gender_cat == "trans" ~ "m",
                              sex == "F" & gender_cat == "cis" ~ "w",
                              sex == "M" & gender_cat == "trans" ~ "w",
                              gender_cat == "non_binary" ~ "nb"),
           trans = if_else(gender_cat == "trans","t","nt"),
           gen_sex_cat = lgbt_cat_shr * gender_cat_shr * pop
           ) %>%
    # sex and acs_bin_col are kept (not just acs_bin) for the cell-grain
    # PUMA calibration below, which needs birth sex (for the singles/married
    # signals) and acs_bin_col (the grain those signals are keyed on).
    group_by(GEOID,puma_id,sex,gender,lgbt_cat,trans,acs_bin,acs_bin_col) %>%
    summarize(total_pop = first(pop),
              gen_sex_cat = sum(gen_sex_cat,na.rm=TRUE))

  # Sum by gender across census block groups
  gender_totals = long2 %>%
    group_by(GEOID,puma_id,gender) %>%
    summarize(gender_pop = sum(gen_sex_cat, na.rm=TRUE)) %>%
    mutate(col_name = paste0("total_",gender)) %>%
    select(-gender) %>%
    pivot_wider(names_from = col_name,
                values_from = gender_pop,
                values_fill = 0)

  # Total adult tract population, computed directly from the pre-split ACS
  # sex x acs_bin population so it's correct regardless of how the
  # cell-level calibration below collapses/regroups acs_bin.
  tract_total_pop <- long %>%
    dplyr::filter(pop > 0) %>%
    dplyr::distinct(GEOID, puma_id, sex, acs_bin, pop) %>%
    dplyr::group_by(GEOID, puma_id) %>%
    dplyr::summarise(total_pop = sum(pop, na.rm = TRUE), .groups = "drop")

  tract_stats <- long2 %>%
    ungroup() %>%
    mutate(
      lgbt_cat = tolower(lgbt_cat),
      col_name = paste(gender, trans, lgbt_cat, sep = "_")
    ) %>%
    # This raw/uncalibrated table doesn't distinguish birth sex within the
    # nb bucket (unlike the cell-grain calibration below, which does) -- sum
    # it away here before pivoting, since long2 now carries sex as its own
    # grouping column and nb cells span both birth sexes.
    group_by(GEOID, puma_id, acs_bin, col_name) %>%
    summarize(total_pop = first(total_pop), gen_sex_cat = sum(gen_sex_cat, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from  = col_name,
      values_from = gen_sex_cat,
      values_fill = 0
    ) %>%
    group_by(GEOID, puma_id) %>%
    summarize(
      total_pop = sum(total_pop, na.rm = TRUE),
      across(matches("^(m|w|nb)_"), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    left_join(gender_totals, by = c("GEOID", "puma_id")) %>%
    mutate(lgbt_m  = m_nt_lg + m_t_lg + m_nt_bisexual + m_t_bisexual + m_nt_queer + m_t_queer + m_t_straight,
           lgbt_w  = w_nt_lg + w_t_lg + w_nt_bisexual + w_t_bisexual + w_nt_queer + w_t_queer + w_t_straight,
           lgbt_nb = nb_nt_lg + nb_nt_bisexual + nb_nt_queer,
           est_lgbt_uncal = lgbt_m + lgbt_w + lgbt_nb,
           lg_m  = m_nt_lg + m_t_lg,
           lg_w  = w_nt_lg + w_t_lg,
           lg_nb = nb_nt_lg,
           bi_m  = m_nt_bisexual + m_t_bisexual,
           bi_w  = w_nt_bisexual + w_t_bisexual,
           bi_nb = nb_nt_bisexual,
           queer_m  = m_nt_queer + m_t_queer,
           queer_w  = w_nt_queer + w_t_queer,
           queer_nb = nb_nt_queer,
           trans_m  = m_t_lg + m_t_bisexual + m_t_queer + m_t_straight,
           trans_w  = w_t_lg + w_t_bisexual + w_t_queer + w_t_straight)
           #est_lgbt_uncal = w_lgbt + m_lgbt + nb_lgbt,
           #est_lg_uncal = w_nt_lg + w_t_lg + m_nt_lg + m_t_lg + nb_nt_lg,
           #est_bi_uncal = w_nt_bisexual + w_t_bisexual + m_nt_bisexual + m_t_bisexual + nb_nt_bisexual,
           #est_queer_uncal = w_nt_queer + w_t_queer + m_nt_queer + m_t_queer + nb_nt_queer,
           #est_trans_uncal = w_t_lg + w_t_bisexual + w_t_queer + w_t_straight + m_t_lg + m_t_bisexual + m_t_queer + m_t_straight)
           #shr_lgbt_uncal = est_lgbt_uncal/total_pop)
  
  #tract_stats <- long2 %>%
  #  group_by(GEOID, puma_id) %>%
  #  summarise(
  #     = sum(pop, na.rm = TRUE),
  #    est_lgbt_uncal = sum(expected_lgbt, na.rm = TRUE),
  #    est_lg_uncal   = sum(expected_lg, na.rm = TRUE),
  #    est_bi_uncal   = sum(expected_bi, na.rm = TRUE),
  #   est_queer_uncal= sum(expected_queer, na.rm = TRUE),
  #    est_trans_uncal= sum(expected_trans, na.rm = TRUE),
  #    .groups = "drop"
  #  ) %>%
  #  mutate(
  #    shr_lgbt_uncal = est_lgbt_uncal / pmax(, 1)
  #  )
  
  # --- Cell-grain PUMA-level calibration ---
  #
  # Generalizes the original single-signal, PUMA-level tilt into a blend of
  # two spatial signals, with the blend weight varying by calibration cell
  # (birth sex x gender bucket x trans status x LGBT category x ACS age
  # bin). Same-sex couples' PUMA concentration is trusted more for cells
  # that are more often actually coupled (per HPS's married-share, the best
  # available proxy); the complementary weight goes to the PUMA
  # concentration of singles of that cell's birth sex/age. See
  # docs/methods.tex for the rationale -- same-sex couples are a biased
  # (older, more legibly partnered) proxy for the full LGBTQ population,
  # worst for exactly the young/non-binary/trans cells least likely to be
  # in a Census-visible couple.
  if (isTRUE(use_calibration)) {

    couples_signal <- get_state_couples_signal(
      state_abbr = state_abbr,
      year = year
    )
    singles_signal <- get_state_singles_signal(
      state_abbr = state_abbr,
      year = year
    )
    married_share <- readRDS(file.path(root_dir, "data/hps", "married_share_by_cell.rds"))

    if (!is.null(couples_signal)) {

      # Fine-grain (tract-level) population per calibration cell. `sex`
      # (birth sex) is redundant with (gender, trans) for m/w buckets but
      # necessary for nb, where AMAB/AFAB isn't otherwise recoverable once
      # gender_cat collapses both into "non_binary".
      tract_cell_pop <- long2 %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          lgbt_cat_lower = tolower(lgbt_cat),
          gen_cat = dplyr::case_when(
            gender == "nb" ~ "non_binary",
            trans == "t"   ~ "trans",
            TRUE            ~ "cis"
          )
        ) %>%
        dplyr::group_by(GEOID, puma_id, sex, gender, trans, lgbt_cat, lgbt_cat_lower, gen_cat, acs_bin_col) %>%
        dplyr::summarise(pop = sum(gen_sex_cat, na.rm = TRUE), .groups = "drop") %>%
        # sex included here (not just gender/trans) so nb's AMAB and AFAB
        # populations calibrate as distinct cells -- redundant with
        # (gender, trans) for m/w but necessary for nb, which otherwise
        # collides on the same cell_id for two different S_uncal/blend
        # values.
        dplyr::mutate(cell_id = paste(sex, gender, trans, lgbt_cat_lower, acs_bin_col, sep = "__"))

      # Same cells, aggregated to PUMA level -- the unit the couples/singles
      # spatial signals and the tilt/gamma calibration operate on.
      puma_cell_pop <- tract_cell_pop %>%
        dplyr::group_by(puma_id, sex, gender, trans, lgbt_cat, lgbt_cat_lower, gen_cat, acs_bin_col) %>%
        dplyr::summarise(pop_uncal = sum(pop, na.rm = TRUE), .groups = "drop") %>%
        dplyr::group_by(sex, gender, trans, lgbt_cat, lgbt_cat_lower, gen_cat, acs_bin_col) %>%
        dplyr::mutate(S_uncal = pop_uncal / sum(pop_uncal, na.rm = TRUE)) %>%
        dplyr::ungroup()

      # Per-cell blend weight: nationwide HPS married-share (the closest
      # available proxy for "coupled", per HPS's MS field -- see
      # code/01_process_hps.R). Floored/ceilinged there to [0.05, 0.95] so
      # no cell relies entirely on one spatial signal.
      puma_cell_pop <- puma_cell_pop %>%
        dplyr::left_join(
          married_share %>% dplyr::select(sex, gen_cat, lgbt_cat, acs_bin_col, married_shr),
          by = c("sex", "gen_cat", "lgbt_cat", "acs_bin_col")
        ) %>%
        dplyr::mutate(w_cell = dplyr::coalesce(married_shr, 0.5))

      # Blended PUMA-level target share per cell. Same-sex-couples
      # concentration is state-level only (not cell-specific, per the
      # existing couples signal); singles concentration varies by this
      # cell's birth sex and age bin. Where only one signal has data for a
      # given PUMA (e.g. filtered out by min_couples/min_singles), fall
      # back to that one signal rather than treating the missing one as a
      # hard zero.
      couple_lookup <- couples_signal %>% dplyr::select(puma_id, couple_share_norm)
      singles_lookup <- if (!is.null(singles_signal)) {
        singles_signal %>% dplyr::select(puma_id, sex, acs_bin_col, singles_share_norm)
      } else {
        tibble::tibble(puma_id = character(), sex = character(),
                        acs_bin_col = character(), singles_share_norm = double())
      }

      puma_cell_pop <- puma_cell_pop %>%
        dplyr::left_join(couple_lookup, by = "puma_id") %>%
        dplyr::left_join(singles_lookup, by = c("puma_id", "sex", "acs_bin_col")) %>%
        dplyr::mutate(
          has_couple  = !is.na(couple_share_norm),
          has_singles = !is.na(singles_share_norm),
          blended_share = dplyr::case_when(
            has_couple & has_singles  ~ w_cell * couple_share_norm + (1 - w_cell) * singles_share_norm,
            has_couple & !has_singles ~ couple_share_norm,
            !has_couple & has_singles ~ singles_share_norm,
            TRUE                      ~ NA_real_
          )
        )

      # Tilt/gamma calibration -- identical math to the original PUMA-level
      # version (gamma keeps its role as a global conservative dial), just
      # computed within each cell instead of once globally. A PUMA with no
      # signal at all for a cell gets blended_share = 0.5 (logit 0, no tilt).
      puma_cell_cal <- puma_cell_pop %>%
        dplyr::mutate(
          blended_share = dplyr::coalesce(blended_share, 0.5),
          blended_share = pmin(pmax(blended_share, 1e-6), 1 - 1e-6),
          puma_signal = qlogis(blended_share)
        ) %>%
        dplyr::group_by(sex, gender, trans, lgbt_cat, lgbt_cat_lower, gen_cat, acs_bin_col) %>%
        dplyr::mutate(
          puma_lgbt_signal = puma_signal - stats::weighted.mean(puma_signal, w = S_uncal),
          tilt_lgbt = exp(gamma * puma_lgbt_signal),
          S_target = (S_uncal * tilt_lgbt) / sum(S_uncal * tilt_lgbt),
          calib_factor = S_target / pmax(S_uncal, 1e-12)
        ) %>%
        dplyr::ungroup() %>%
        # sex included here (not just gender/trans) so nb's AMAB and AFAB
        # populations calibrate as distinct cells -- redundant with
        # (gender, trans) for m/w but necessary for nb, which otherwise
        # collides on the same cell_id for two different S_uncal/blend
        # values.
        dplyr::mutate(cell_id = paste(sex, gender, trans, lgbt_cat_lower, acs_bin_col, sep = "__")) %>%
        dplyr::select(puma_id, cell_id, calib_factor)

      # IDW-smooth calibration factors to GEOID level so hard PUMA boundary
      # discontinuities are replaced by smooth spatial transitions. One
      # distance/weight matrix is shared across every cell -- geometry
      # doesn't depend on cell -- applied via a single matrix multiply
      # against a (PUMA x cell) factor matrix.
      #
      # Centroids and distances are computed in geographic coordinates
      # (EPSG:4326) rather than projected Web Mercator (EPSG:3857). Mercator
      # has no antimeridian wraparound, so a naive projected distance between
      # a tract at e.g. 179.5 E and one at 179.5 W treats them as being on
      # opposite sides of the world instead of ~1 degree apart -- this
      # silently corrupted IDW weights for Alaska's western Aleutian Islands
      # tracts. With sf_use_s2() enabled (the default here), st_centroid()
      # and st_distance() on geographic geometries use spherical (S2)
      # geometry, which handles the antimeridian correctly and also avoids
      # Mercator's latitude-dependent distance distortion (significant this
      # far north regardless of the antimeridian issue).
      puma_cell_wide <- puma_cell_cal %>%
        tidyr::pivot_wider(names_from = cell_id, values_from = calib_factor, values_fill = 1)

      cell_cols <- setdiff(names(puma_cell_wide), "puma_id")

      puma_cents <- pumas_clean |>
        dplyr::inner_join(puma_cell_wide, by = "puma_id") |>
        sf::st_transform(4326) |>
        sf::st_centroid()

      tract_cents <- tract_clean |>
        dplyr::filter(GEOID %in% unique(tract_cell_pop$GEOID)) |>
        dplyr::select(GEOID) |>
        sf::st_transform(4326) |>
        sf::st_centroid()

      dist_mat <- sf::st_distance(tract_cents, puma_cents)
      dist_num <- matrix(as.numeric(dist_mat), nrow = nrow(tract_cents))
      idw_w    <- 1 / pmax(dist_num, 1)^2
      idw_w    <- idw_w / rowSums(idw_w)

      puma_cell_matrix <- as.matrix(sf::st_drop_geometry(puma_cents)[, cell_cols, drop = FALSE])
      tract_cell_factor <- idw_w %*% puma_cell_matrix
      colnames(tract_cell_factor) <- cell_cols

      tract_cell_idw <- tibble::as_tibble(tract_cell_factor) %>%
        dplyr::mutate(GEOID = tract_cents$GEOID) %>%
        tidyr::pivot_longer(cols = dplyr::all_of(cell_cols), names_to = "cell_id", values_to = "calib_factor")

      # Apply calibration at cell grain, then collapse (sex, acs_bin_col)
      # back into the gender/trans/lgbt_cat subcategory columns the rest of
      # the pipeline (and app.R) expects -- each subcategory is now built
      # from its own per-cell-calibrated population rather than one flat
      # PUMA factor applied uniformly.
      tract_stats_cal <- tract_cell_pop %>%
        dplyr::left_join(tract_cell_idw, by = c("GEOID", "cell_id")) %>%
        dplyr::mutate(
          calib_factor = dplyr::coalesce(calib_factor, 1),
          pop_map = pop * calib_factor
        ) %>%
        dplyr::group_by(GEOID, puma_id, gender, trans, lgbt_cat_lower) %>%
        dplyr::summarise(pop_map = sum(pop_map, na.rm = TRUE), .groups = "drop") %>%
        dplyr::mutate(col_name = paste(gender, trans, lgbt_cat_lower, sep = "_")) %>%
        dplyr::select(GEOID, puma_id, col_name, pop_map) %>%
        tidyr::pivot_wider(names_from = col_name, values_from = pop_map, values_fill = 0) %>%
        dplyr::left_join(tract_total_pop, by = c("GEOID", "puma_id")) %>%
        dplyr::left_join(gender_totals, by = c("GEOID", "puma_id")) %>%
        dplyr::mutate(
          lgbt_m_map_raw   = m_nt_lg + m_t_lg + m_nt_bisexual + m_t_bisexual + m_nt_queer + m_t_queer + m_t_straight,
          lgbt_w_map_raw   = w_nt_lg + w_t_lg + w_nt_bisexual + w_t_bisexual + w_nt_queer + w_t_queer + w_t_straight,
          lgbt_nb_map_raw  = nb_nt_lg + nb_nt_bisexual + nb_nt_queer,
          lg_m_map_raw     = m_nt_lg + m_t_lg,
          lg_w_map_raw     = w_nt_lg + w_t_lg,
          lg_nb_map_raw    = nb_nt_lg,
          bi_m_map_raw     = m_nt_bisexual + m_t_bisexual,
          bi_w_map_raw     = w_nt_bisexual + w_t_bisexual,
          bi_nb_map_raw    = nb_nt_bisexual,
          queer_m_map_raw  = m_nt_queer + m_t_queer,
          queer_w_map_raw  = w_nt_queer + w_t_queer,
          queer_nb_map_raw = nb_nt_queer,
          trans_m_map_raw  = m_t_lg + m_t_bisexual + m_t_queer + m_t_straight,
          trans_w_map_raw  = w_t_lg + w_t_bisexual + w_t_queer + w_t_straight
        )

      scale_lgbt <-
        sum(tract_stats$lgbt_m + tract_stats$lgbt_w + tract_stats$lgbt_nb, na.rm = TRUE) /
        sum(tract_stats_cal$lgbt_m_map_raw + tract_stats_cal$lgbt_w_map_raw + tract_stats_cal$lgbt_nb_map_raw, na.rm = TRUE)

      tract_stats_cal <- tract_stats_cal %>%
        mutate(
          lgbt_m_map   = lgbt_m_map_raw  * scale_lgbt,
          lgbt_w_map   = lgbt_w_map_raw  * scale_lgbt,
          lgbt_nb_map  = lgbt_nb_map_raw * scale_lgbt,
          lg_m_map     = lg_m_map_raw * scale_lgbt,
          lg_w_map     = lg_w_map_raw * scale_lgbt,
          lg_nb_map    = lg_nb_map_raw * scale_lgbt,
          bi_m_map     = bi_m_map_raw * scale_lgbt,
          bi_w_map     = bi_w_map_raw * scale_lgbt,
          bi_nb_map    = bi_nb_map_raw * scale_lgbt,
          queer_m_map  = queer_m_map_raw * scale_lgbt,
          queer_w_map  = queer_w_map_raw * scale_lgbt,
          queer_nb_map = queer_nb_map_raw * scale_lgbt,
          trans_m_map  = trans_m_map_raw * scale_lgbt,
          trans_w_map  = trans_w_map_raw * scale_lgbt
          #est_lg_map    = est_lg_map_raw    * scale_lgbt,
          #est_bi_map    = est_bi_map_raw    * scale_lgbt,
          #est_queer_map = est_queer_map_raw * scale_lgbt,
          #est_trans_map = est_trans_map_raw * scale_lgbt,
          #shr_lgbt_map  = est_lgbt_map / pmax(total_pop, 1)
        ) |>
        # --- Enforce CBG capacity constraint: est_lgb_map <= pop18 ---
        #dplyr::group_by(puma_id) |>
        #dplyr::mutate(
          # capacity at each CBG
        #  capacity = pmax(total_pop - est_lgbt_map, 0),
          
          # excess mass at each CBG
        #  excess = pmax(est_lgbt_map - total_pop, 0)
        #) |>
        #dplyr::mutate(
        #  total_excess   = sum(excess),
        #  total_capacity = sum(capacity)
        #) |>
        #dplyr::mutate(
          # redistribute excess proportionally to capacity
          #est_lgbt_map = dplyr::if_else(
          #  total_excess > 0 & total_capacity > 0,
          #  pmin(
          #    total_pop,
          #    est_lgbt_map - excess + capacity * (total_excess / total_capacity)
          # ),
          #  est_lgbt_map
          #)
        #) |>
        dplyr::ungroup() |>
        # Cap each gender group so LGBTQ estimates never exceed total population.
        # scale proportionally so within-group composition is preserved.
        dplyr::mutate(
          cap_m  = pmin(1, total_m  / pmax(lgbt_m_map,  1e-9)),
          cap_w  = pmin(1, total_w  / pmax(lgbt_w_map,  1e-9)),
          cap_nb = pmin(1, total_nb / pmax(lgbt_nb_map, 1e-9)),
          lgbt_m_map   = lgbt_m_map   * cap_m,
          lgbt_w_map   = lgbt_w_map   * cap_w,
          lgbt_nb_map  = lgbt_nb_map  * cap_nb,
          lg_m_map     = lg_m_map     * cap_m,
          lg_w_map     = lg_w_map     * cap_w,
          lg_nb_map    = lg_nb_map    * cap_nb,
          bi_m_map     = bi_m_map     * cap_m,
          bi_w_map     = bi_w_map     * cap_w,
          bi_nb_map    = bi_nb_map    * cap_nb,
          queer_m_map  = queer_m_map  * cap_m,
          queer_w_map  = queer_w_map  * cap_w,
          queer_nb_map = queer_nb_map * cap_nb,
          trans_m_map  = trans_m_map  * cap_m,
          trans_w_map  = trans_w_map  * cap_w
        ) |>
        dplyr::mutate(
          shr_lgbt_m_map = lgbt_m_map / pmax(total_m, 1),
          shr_lgbt_w_map = lgbt_w_map / pmax(total_w, 1),
          shr_lgbt_nb_map = lgbt_nb_map / pmax(total_nb, 1),
          shr_lg_m_map = lg_m_map / pmax(total_m, 1),
          shr_lg_w_map = lg_w_map / pmax(total_w, 1),
          shr_lg_nb_map = lg_nb_map / pmax(total_nb, 1),
          shr_bi_m_map = bi_m_map / pmax(total_m, 1),
          shr_bi_w_map = bi_w_map / pmax(total_w, 1),
          shr_bi_nb_map = bi_nb_map / pmax(total_nb, 1),
          shr_queer_m_map = queer_m_map / pmax(total_m, 1),
          shr_queer_w_map = queer_w_map / pmax(total_w, 1),
          shr_queer_nb_map = queer_nb_map / pmax(total_nb, 1),
          shr_trans_m_map = trans_m_map / pmax(total_m, 1),
          shr_trans_w_map = trans_w_map / pmax(total_w, 1),
        )
      
    } 
  }
  
  tract_puma_xwalk |>
    dplyr::left_join(tract_stats_cal, by = c("GEOID","puma_id")) |>
    sf::st_transform(CRS_LEAFLET)
}
