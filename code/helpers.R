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
      puma_signal = qlogis(p)
    ) |>
    dplyr::select(puma_id, puma_signal)
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

    fallback_match <- sf::st_join(
      fallback_geoms, cbsas_sf, join = sf::st_intersects, largest = TRUE
    ) |>
      sf::st_drop_geometry() |>
      distinct(GEOID, .keep_all = TRUE) |>
      rename(cbsa_id_fallback = cbsa_id, cbsa_name_fallback = cbsa_name)

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
    group_by(GEOID,puma_id,gender,lgbt_cat,trans,acs_bin) %>%
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
  
  tract_stats <- long2 %>%
    ungroup() %>%
    mutate(
      lgbt_cat = tolower(lgbt_cat),
      col_name = paste(gender, trans, lgbt_cat, sep = "_")
    ) %>%
    select(GEOID, puma_id, acs_bin, total_pop, col_name, gen_sex_cat) %>%
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
  
  # --- PUMA-level calibration ---
  if (isTRUE(use_calibration)) {
    
    signal <- get_state_couples_signal(
      state_abbr = state_abbr,
      year = year
    )
    
    if (!is.null(signal)) {
      
      puma_uncal <- tract_stats %>%
        group_by(puma_id) %>%
        summarise(
          lgbt_uncal = sum(est_lgbt_uncal, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(S_lgbt_uncal = lgbt_uncal / sum(lgbt_uncal, na.rm = TRUE))
      
      puma_tilt <- puma_uncal |>
        dplyr::left_join(signal, by = "puma_id") |>
        dplyr::mutate(
          puma_signal = dplyr::coalesce(puma_signal, 0),
          puma_lgbt_signal = puma_signal - stats::weighted.mean(puma_signal, w = S_lgbt_uncal),
          tilt_lgbt = exp(gamma * puma_lgbt_signal),
          S_lgbt_target = (S_lgbt_uncal * tilt_lgbt) / sum(S_lgbt_uncal * tilt_lgbt),
          calib_lgbt_factor = S_lgbt_target / pmax(S_lgbt_uncal, 1e-12)
        ) |>
        dplyr::select(puma_id, calib_lgbt_factor) %>%
        distinct(puma_id, .keep_all = TRUE)
      
      # IDW-smooth calibration factors to GEOID level so hard PUMA boundary
      # discontinuities are replaced by smooth spatial transitions.
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
      puma_cents <- pumas_clean |>
        dplyr::inner_join(
          dplyr::select(puma_tilt, puma_id, calib_lgbt_factor),
          by = "puma_id"
        ) |>
        sf::st_transform(4326) |>
        sf::st_centroid()

      tract_cents <- tract_clean |>
        dplyr::filter(GEOID %in% unique(tract_stats$GEOID)) |>
        dplyr::select(GEOID) |>
        sf::st_transform(4326) |>
        sf::st_centroid()

      dist_mat <- sf::st_distance(tract_cents, puma_cents)
      dist_num <- matrix(as.numeric(dist_mat), nrow = nrow(tract_cents))
      idw_w    <- 1 / pmax(dist_num, 1)^2
      idw_w    <- idw_w / rowSums(idw_w)

      tract_idw <- tibble::tibble(
        GEOID             = tract_cents$GEOID,
        calib_lgbt_factor = as.numeric(idw_w %*% puma_cents$calib_lgbt_factor)
      )

      tract_stats_cal <- tract_stats |>
        dplyr::left_join(tract_idw, by = "GEOID") |>
        dplyr::mutate(
          calib_lgbt_factor = dplyr::coalesce(calib_lgbt_factor, 1),
          lgbt_m_map_raw   = lgbt_m  * calib_lgbt_factor,
          lgbt_w_map_raw   = lgbt_w  * calib_lgbt_factor,
          lgbt_nb_map_raw  = lgbt_nb * calib_lgbt_factor,
          lg_m_map_raw     = lg_m * calib_lgbt_factor,
          lg_w_map_raw     = lg_w * calib_lgbt_factor,
          lg_nb_map_raw    = lg_nb * calib_lgbt_factor,
          bi_m_map_raw     = bi_m * calib_lgbt_factor,
          bi_w_map_raw     = bi_w * calib_lgbt_factor,
          bi_nb_map_raw   = bi_nb * calib_lgbt_factor,
          queer_m_map_raw  = queer_m * calib_lgbt_factor,
          queer_w_map_raw  = queer_w * calib_lgbt_factor,
          queer_nb_map_raw = queer_nb * calib_lgbt_factor,
          trans_m_map_raw  = trans_m * calib_lgbt_factor,
          trans_w_map_raw  = trans_w * calib_lgbt_factor
          #est_lgbt_map_raw  = est_lgbt_uncal  * calib_lgbt_factor,
          #est_lg_map_raw    = est_lg_uncal    * calib_lgbt_factor,
          #est_bi_map_raw    = est_bi_uncal    * calib_lgbt_factor,
          #est_queer_map_raw = est_queer_uncal * calib_lgbt_factor,
          #est_trans_map_raw = est_trans_uncal * calib_lgbt_factor
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
