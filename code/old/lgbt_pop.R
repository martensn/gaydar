library(tidycensus)
library(tigris)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(units)
library(tidygeocoder)

options(tigris_use_cache = TRUE)
CRS_US <- 5070  # NAD83 / Conus Albers (meters)


file_list <- list.files("/Users/nickmartens/Desktop/Goldenberg/LGBTQ/Raw/")
file_list <- file_list[str_detect(file_list, ".csv") & !str_detect(file_list, "repwgt")]

# bind together
hps <- paste0("/Users/nickmartens/Desktop/Goldenberg/LGBTQ/Raw/", file_list) %>% 
  lapply(fread) %>% 
  bind_rows()

ref_year <- 2023
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

hps_acs_rates = hps %>%
  as_tibble() %>%
  mutate(sexuality = case_when(SEXUAL_ORIENTATION == 2 ~ "Straight",
                               SEXUAL_ORIENTATION == 3 ~ "Bisexual",
                               SEXUAL_ORIENTATION %in% c(1,4:5) ~ "LG",
                               TRUE ~ "Straight"),
         sex = case_when(EGENID_BIRTH == 1 ~ "M",
                         EGENID_BIRTH == 2 ~ "F"),
         age = ref_year - TBIRTH_YEAR) %>%
  #SEXUAL_ORIENTATION == 3 ~ "Bisexual",
  #SEXUAL_ORIENTATION == 4 ~ "Other Queer",
  #SEXUAL_ORIENTATION == 5 ~ "Questioning"),
  #Gender = case_when(GENID_DESCRIBE == 1 ~ "AMAB",
  #GENID_DESCRIBE == 2 ~ "AFAB",
  #GENID_DESCRIBE == 3 ~ "Trans",
  #GENID_DESCRIBE == 4 ~ "GNC")) %>%
  group_by(sex,sexuality,EST_ST,age,WEEK) %>%
  summarize(total = sum(PWEIGHT)) %>%
  left_join(
    acs_age_bins,
    by = join_by(age >= age_lo, age <= age_hi)
  ) %>%
  group_by(sex,sexuality,EST_ST,acs_bin,WEEK) %>%
  summarize(total = sum(total)) %>%
  group_by(sex,sexuality,EST_ST,acs_bin) %>%
  summarize(total = mean(total,na.rm=TRUE)) %>%
  pivot_wider(names_from=sexuality,
              values_from=total,
              values_fill = 0) %>%
  mutate(lg_shr = LG/(LG+Straight+Bisexual),
         lgb_shr = (LG+Bisexual)/(LG+Straight+Bisexual)) %>%
  left_join(state_crosswalk, by = c("EST_ST"="state_fips")) %>%
  select(sex,state,state_abbr,acs_bin,lg_shr,lgb_shr)





get_bg_in_buffer <- function(buf, year = 2023) {
  
  # identify intersecting counties dynamically
  counties_hit <- tigris::counties(cb = TRUE, year = year) |>
    st_transform(CRS_US) |>
    st_intersects(buf, sparse = FALSE) |>
    (\(m) dplyr::filter(
      tigris::counties(cb = TRUE, year = year),
      m[,1]
    ))()
  
  # pull B01001 for those counties
  pmap_dfr(
    list(counties_hit$STATEFP, counties_hit$COUNTYFP),
    \(statefp, countyfp) {
      get_acs(
        geography = "block group",
        variables = c(
          # total pop isn’t strictly required if you sum age-sex bins, but handy
          total = "B01001_001",
          # MALES
          m_0_4  = "B01001_003", m_5_9  = "B01001_004", m_10_14="B01001_005",
          m_15_17="B01001_006", m_18_19="B01001_007", m_20  ="B01001_008",
          m_21  ="B01001_009", m_22_24="B01001_010", m_25_29="B01001_011",
          m_30_34="B01001_012", m_35_39="B01001_013", m_40_44="B01001_014",
          m_45_49="B01001_015", m_50_54="B01001_016", m_55_59="B01001_017",
          m_60_61="B01001_018", m_62_64="B01001_019", m_65_66="B01001_020",
          m_67_69="B01001_021", m_70_74="B01001_022", m_75_79="B01001_023",
          m_80_84="B01001_024", m_85up ="B01001_025",
          # FEMALES
          f_0_4  ="B01001_027", f_5_9  ="B01001_028", f_10_14="B01001_029",
          f_15_17="B01001_030", f_18_19="B01001_031", f_20  ="B01001_032",
          f_21  ="B01001_033", f_22_24="B01001_034", f_25_29="B01001_035",
          f_30_34="B01001_036", f_35_39="B01001_037", f_40_44="B01001_038",
          f_45_49="B01001_039", f_50_54="B01001_040", f_55_59="B01001_041",
          f_60_61="B01001_042", f_62_64="B01001_043", f_65_66="B01001_044",
          f_67_69="B01001_045", f_70_74="B01001_046", f_75_79="B01001_047",
          f_80_84="B01001_048", f_85up ="B01001_049"
        ),   # define once globally
        state = statefp,
        county = countyfp,
        year = year,
        geometry = TRUE,
        output = "wide",
        survey = "acs5"
      )
    }
  ) |>
    st_transform(CRS_US)
}

area_weight_pop <- function(bgs, buf) {
  
  bgs <- bgs |> mutate(area_bg = st_area(geometry))
  
  bgs_in <- st_intersection(bgs, buf) |>
    mutate(
      area_in = st_area(geometry),
      w = as.numeric(area_in / area_bg)
    )
  
  bgs_in |>
    st_drop_geometry() |>
    transmute(
      across(where(is.numeric) & ends_with("E"), ~ .x * w)
    ) |>
    summarise(across(where(is.numeric), sum, na.rm = TRUE))
}

acs_to_long <- function(pop_df) {
  pop_df |>
    pivot_longer(
      cols = matches("^[mf]_"),
      names_to = "var",
      values_to = "pop"
    ) |>
    mutate(
      sex = if_else(str_starts(var, "m_"), "M", "F"),
      acs_bin = var |>
        str_remove("^[mf]_") |>
        str_remove("E$")
    ) |>
    select(sex, acs_bin, pop)
}

estimate_lgb_for_address <- function(address_row,
                                     radius_miles = 10,
                                     year = 2023) {
  
  # 1) geocode
  pt <- tidygeocoder::geocode(
    address_row,
    street = street,
    city = city,
    state = state,
    postalcode = zip,
    method = "census"
  ) |>
    st_as_sf(coords = c("long","lat"), crs = 4326) |>
    st_transform(CRS_US)
  
  # 2) buffer
  buf <- st_buffer(pt, dist = units::set_units(radius_miles, "miles"))
  
  # 3) pull BGs
  bgs <- get_bg_in_buffer(buf, year)
  
  # 4) area-weighted population
  pop_in_radius <- area_weight_pop(bgs, buf)
  
  # 5) long ACS
  acs_long <- acs_to_long(pop_in_radius)
  
  # -------------------------------
  # 6A) PURE ACS POPULATION STATS
  # -------------------------------
  pop_stats <- acs_long |>
    mutate(
      is_22_29 = acs_bin %in% c("22_24", "25_29"),
      is_male  = sex == "M",
      pop_22_29   = case_when(is_22_29 == TRUE ~ pop,
                              TRUE ~ 0),
      pop_m_22_29 = case_when(is_22_29 == TRUE & is_male == TRUE ~ pop,
                              TRUE ~ 0)
    ) |>
    summarise(
      pop = sum(pop, na.rm = TRUE),
      pop_22_29 = sum(pop_22_29, na.rm=TRUE),
      pop_m_22_29 = sum(pop_m_22_29, na.rm=TRUE)
    )
  
  # --------------------------------
  # 6B) EXPECTED LG / LGB COUNTS
  # --------------------------------
  rates_state <- hps_acs_rates |>
    filter(state_abbr == address_row$state) |>
    group_by(sex, acs_bin) |>
    summarise(
      lg_shr  = mean(lg_shr,  na.rm = TRUE),
      lgb_shr = mean(lgb_shr, na.rm = TRUE),
      .groups = "drop"
    )
  
  est_stats <- acs_long |>
    left_join(rates_state, by = c("sex", "acs_bin")) |>
    mutate(
      expected_lg  = pop * lg_shr,
      expected_lgb = pop * lgb_shr,
      is_22_29 = acs_bin %in% c("22_24", "25_29"),
      is_male  = sex == "M"
    ) |>
    summarise(
      est_lg  = sum(expected_lg,  na.rm = TRUE),
      est_lgb = sum(expected_lgb, na.rm = TRUE),
      
      est_lg_22_29   = sum(if_else(is_22_29, expected_lg,  0), na.rm = TRUE),
      est_lgb_22_29  = sum(if_else(is_22_29, expected_lgb, 0), na.rm = TRUE),
      
      est_lg_m_22_29 = sum(if_else(is_22_29 & is_male, expected_lg,  0), na.rm = TRUE),
      est_lg_f_22_29 = sum(if_else(is_22_29 & !is_male, expected_lg, 0), na.rm = TRUE),
      est_lgb_m_22_29 = sum(if_else(is_22_29 & is_male, expected_lgb, 0), na.rm = TRUE),
      
      .groups = "drop"
    )
  
  # -------------------------------
  # 6C) COMBINE + FINAL METRICS
  # -------------------------------
  dplyr::bind_cols(pop_stats, est_stats) |>
    mutate(
      shr_fruity = est_lgb_m_22_29 / pop_m_22_29
    )
}


# tidycensus API key (do once)
# census_api_key("YOUR_KEY_HERE", install = TRUE)
address <- tibble::tribble(
  ~street, ~city, ~state, ~zip, 
  "1645 E 50th St", "Chicago", "IL", "60615",
  "380 Thistle Ave", "Climax", "MI", "49034",
  "259 Nassau St", "Princeton", "NJ", "08542",
  "611 Tappan Ave", "Ann Arbor", "MI", "48109"
)

lgbt_results <- address |>
  mutate(
    stats = purrr::map(
      seq_len(n()),
      ~ estimate_lgb_for_address(address[.x, ])
    )
  ) |>
  tidyr::unnest(stats)
