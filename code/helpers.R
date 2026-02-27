# R/helpers.R
# Utilities for building a BG-level sf layer with expected LGB counts/shares.
library(sf)
library(dplyr)
library(tidyr)
library(stringr)
library(tidycensus)

CRS_LEAFLET <- 4326

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

get_state_choices <- function() {
  tigris::states(cb = TRUE, year = 2023) |>
    sf::st_drop_geometry() |>
    dplyr::arrange(NAME) |>
    dplyr::pull(STUSPS)
}



# Convert ACS wide B01001 into long with sex + acs_bin for adults only
acs_wide_to_long_adults <- function(bg_wide, bg_puma_xwalk) {
  # bg_wide: sf with columns matching names(B01001_VARS), including GEOID and geometry
  long <- bg_wide |>
    sf::st_drop_geometry() |>
    inner_join(bg_puma_xwalk, by = "GEOID", relationship = "one-to-many") %>%
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
    dplyr::select(GEOID, puma_id, sex, acs_bin, pop, w_puma)

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


build_bg_expected_layer <- function(
    state_abbr,
    year,
    rates,
    use_calibration = TRUE,
    gamma = 1
) {
  
  # --- BG-level age-sex data ---
  acs_raw <- tidycensus::get_acs(
    geography = "block group",
    variables = B01001_VARS,
    state = state_abbr,
    year = year,
    survey = "acs5",
    geometry = TRUE,
    output = "tidy",
    cache_table = TRUE
  )
  
  bg_wide <- acs_raw |>
    dplyr::select(GEOID, variable, estimate, geometry) |>
    tidyr::pivot_wider(
      names_from = "variable",
      values_from = "estimate",
      values_fill = 0
    ) |>
    sf::st_as_sf()
  
  # --- Map BG → PUMA (spatial join once) ---
  pumas <- tigris::pumas(
    state = state_abbr,
    year = 2019,
    cb = TRUE
  ) |>
    dplyr::select(puma_id = GEOID10)
  
  cbsas <- tigris::core_based_statistical_areas(
    year = 2019,
    cb = TRUE
  ) |>
    dplyr::select(cbsa_id = GEOID, cbsa_name = NAME) |>
    st_drop_geometry()
  
  county_cbsa_xwalk = read_csv(file.path(root_dir,"data/unified_cbsa.csv")) %>%
    select(cbsa_id = cbsa_code, county_fips = GeoFIPS) %>%
    left_join(cbsas, by = "cbsa_id")
  
  bg_clean <- bg_wide |>
    st_transform(st_crs(pumas)) |>
    sf::st_make_valid()
  
  pumas_clean <- pumas |>
    sf::st_make_valid()
  
  bg_puma_xwalk <- st_intersection(
    bg_clean |> select(GEOID),
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
  
  # --- Adult population long ---
  long <- acs_wide_to_long_adults(bg_wide, bg_puma_xwalk)
  
  rates <- readRDS(file.path(root_dir,"data/hps/hps_acs_rates.rds"))
  
  rates_state <- rates %>%
    filter(state_abbr == !!state_abbr) %>%
    select(-state_abbr)
  
  long2 <- long %>%
    left_join(rates_state, by = c("sex", "acs_bin")) %>%
    mutate(
      expected_lgbt = pop * lgbt_shr,
      expected_lg   = pop * coalesce(LG, 0),
      expected_bi   = pop * coalesce(Bisexual, 0),
      expected_queer= pop * coalesce(Queer, 0),
      expected_trans= pop * coalesce(Trans, 0)
    )
  
  bg_stats <- long2 %>%
    group_by(GEOID, puma_id) %>%
    summarise(
      pop18 = sum(pop, na.rm = TRUE),
      est_lgbt_uncal = sum(expected_lgbt, na.rm = TRUE),
      est_lg_uncal   = sum(expected_lg, na.rm = TRUE),
      est_bi_uncal   = sum(expected_bi, na.rm = TRUE),
      est_queer_uncal= sum(expected_queer, na.rm = TRUE),
      est_trans_uncal= sum(expected_trans, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      shr_lgbt_uncal = est_lgbt_uncal / pmax(pop18, 1)
    )
  
  # --- PUMA-level calibration ---
  if (isTRUE(use_calibration)) {
    
    signal <- get_state_couples_signal(
      state_abbr = state_abbr,
      year = year
    )
    
    if (!is.null(signal)) {
      
      puma_uncal <- bg_stats %>%
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
      
      bg_stats2 <- bg_stats %>%
        left_join(puma_tilt %>% select(puma_id, calib_lgbt_factor), by="puma_id") %>%
        mutate(
          calib_lgbt_factor = coalesce(calib_lgbt_factor, 1),
          
          est_lgbt_map_raw  = est_lgbt_uncal  * calib_lgbt_factor,
          est_lg_map_raw    = est_lg_uncal    * calib_lgbt_factor,
          est_bi_map_raw    = est_bi_uncal    * calib_lgbt_factor,
          est_queer_map_raw = est_queer_uncal * calib_lgbt_factor,
          est_trans_map_raw = est_trans_uncal * calib_lgbt_factor
        )
      
      scale_lgbt <- sum(bg_stats2$est_lgbt_uncal, na.rm=TRUE) / sum(bg_stats2$est_lgbt_map_raw, na.rm=TRUE)
      
      bg_stats2 <- bg_stats2 %>%
        mutate(
          est_lgbt_map  = est_lgbt_map_raw  * scale_lgbt,
          est_lg_map    = est_lg_map_raw    * scale_lgbt,
          est_bi_map    = est_bi_map_raw    * scale_lgbt,
          est_queer_map = est_queer_map_raw * scale_lgbt,
          est_trans_map = est_trans_map_raw * scale_lgbt,
          shr_lgbt_map  = est_lgbt_map / pmax(pop18, 1)
        ) |>
        # --- Enforce CBG capacity constraint: est_lgb_map <= pop18 ---
        dplyr::group_by(puma_id) |>
        dplyr::mutate(
          # capacity at each CBG
          capacity = pmax(pop18 - est_lgbt_map, 0),
          
          # excess mass at each CBG
          excess = pmax(est_lgbt_map - pop18, 0)
        ) |>
        dplyr::mutate(
          total_excess   = sum(excess),
          total_capacity = sum(capacity)
        ) |>
        dplyr::mutate(
          # redistribute excess proportionally to capacity
          est_lgbt_map = dplyr::if_else(
            total_excess > 0 & total_capacity > 0,
            pmin(
              pop18,
              est_lgbt_map - excess + capacity * (total_excess / total_capacity)
            ),
            est_lgbt_map
          )
        ) |>
        dplyr::ungroup() |>
        dplyr::mutate(
          shr_lgbt_map = est_lgbt_map / pmax(pop18, 1),
          shr_lg_map = est_lg_map / pmax(pop18, 1),
          shr_bi_map = est_bi_map / pmax(pop18, 1),
          shr_queer_map = est_queer_map / pmax(pop18, 1),
          shr_trans_map = est_trans_map / pmax(pop18, 1)
        )
      
    } 
  }
  
  bg_puma_xwalk |>
    dplyr::left_join(bg_stats2, by = c("GEOID","puma_id")) |>
    sf::st_transform(CRS_LEAFLET)
}


# Identify the cbsa



