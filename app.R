# app.R — Gaydar (address-based comparison dashboard)
# ---------------------------------------------------
library(readr)
library(shiny)
library(sf)
library(dplyr)
library(tidyr)
library(leaflet)
library(ggplot2)
library(tidygeocoder)
library(tidycensus)
library(tigris)
library(htmltools)
library(plotly)
library(rmapshaper)
library(h3jsr)
library(shinycssloaders)

options(tigris_use_cache = TRUE)

source("code/helpers.R")

# ---- national tract-level lookup table (loaded once per process) ----
# Built offline by code/build_national_tract_rates.R from the same cached
# per-state tract layers used elsewhere in the app; used for national
# percentile comparisons and as the source of each tract's Congressional
# District assignment.
national_tract_rates <- readRDS(file.path(root_dir, "data/cache/national_tract_rates.rds"))

# ---- human-readable names for county/PUMA/CD GEOIDs (built offline by
# code/build_geo_names.R) so the map can show e.g. "Alger County, MI"
# instead of a raw FIPS code. Tracts have no such lookup (names = NULL). ----
geo_names <- readRDS(file.path(root_dir, "data/cache/geo_names.rds"))

# ---- geography levels available in the "Geography" selector ----
# col = the column identifying that unit on a tract row; label = display
# name (used in the tooltip/popup); names = key into geo_names, or NULL.
geo_level_info <- list(
  tract  = list(col = "GEOID",       label = "census tract",           names = NULL),
  county = list(col = "county_fips", label = "county",                 names = "county"),
  puma   = list(col = "puma_id",     label = "puma",                   names = "puma"),
  cd     = list(col = "cd_id",       label = "congressional district", names = "cd")
)

# Dissolves a per-state tract sf layer up to county/PUMA/CD by summing the
# numeric count columns within each group (geometry auto-dissolves via sf's
# dplyr method). cbsa_id is kept via first() since it's categorical, needed
# by metro_bg() at every geography level.
dissolve_geo <- function(sf_obj, group_col) {
  sf_obj |>
    dplyr::group_by(.data[[group_col]]) |>
    dplyr::summarise(
      dplyr::across(dplyr::where(is.numeric), \(x) sum(x, na.rm = TRUE)),
      cbsa_id = dplyr::first(cbsa_id),
      .groups = "drop"
    ) |>
    dplyr::rename(GEOID = dplyr::all_of(group_col)) |>
    dplyr::filter(!is.na(GEOID))
}

# Same idea, applied to the geometry-free national lookup table, for
# building the national comparison distribution at the selected geography.
aggregate_national <- function(group_col) {
  if (group_col == "GEOID") return(national_tract_rates)
  national_tract_rates |>
    dplyr::filter(!is.na(.data[[group_col]])) |>
    dplyr::group_by(.data[[group_col]]) |>
    dplyr::summarise(dplyr::across(dplyr::where(is.numeric), \(x) sum(x, na.rm = TRUE)), .groups = "drop")
}

# ---------------------------
# UI
# ---------------------------

ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", href = "app-theme.css"),
    tags$style(HTML("
      html, body {
        height: 100%;
        margin: 0;
      }

      .container-fluid {
        height: 100%;
      }

      .fill-height {
        height: calc(100vh - 80px); /* subtract title height */
        display: flex;
        flex-direction: column;
      }

      .map-panel {
        flex: 1 1 auto;
      }

      .donut-panel {
        flex: 0 0 300px;
      }

      /* shinycssloaders wraps each output in a container div that doesn't
         inherit height, which breaks height:100% on leafletOutput/plotlyOutput
         inside this flex layout (Leaflet in particular renders a blank map if
         its container is 0-height at init) */
      .shiny-spinner-output-container {
        height: 100%;
      }

      @media (max-width: 767px) {
        .fill-height {
          height: auto;
        }
        .map-panel {
          min-height: 400px;
        }
        .donut-panel {
          flex: 0 0 auto;
          min-height: 300px;
        }
      }
    "))
  ),
  titlePanel(title = "Gaydar", windowTitle = "Gaydar — LGBTQ Population Estimator"),
  sidebarLayout(
    sidebarPanel(
      textInput("street", "address", placeholder = "123 Main St"),
      textInput("city",   "city",    placeholder = "Chicago"),
      textInput("state",  "state",   placeholder = "IL"),
      textInput("zip",    "zip",     placeholder = "60637"),
      sliderInput(
        "radius_miles",
        "Radius (miles)",
        min   = 0,
        max   = 50,
        value = 10,
        step  = 1,
        ticks = FALSE
      ),
      selectInput(
        "geo_level",
        "geography",
        choices = c(
          "Census Tract"           = "tract",
          "County"                 = "county",
          "PUMA"                   = "puma",
          "Congressional District" = "cd"
        ),
        selected = "tract"
      ),
      selectInput(
        "metric",
        "population",
        choices = c(
          "LGBT"            = "lgbt",
          "Lesbian and Gay" = "lg",
          "Bisexual"        = "bi",
          "Queer"           = "queer",
          "Transgender"     = "trans"
        ),
        selected = "lgbt"
      ),

      checkboxGroupInput(
        "gender",
        "gender",
        choices = c("Man" = "m", "Woman" = "w", "Non-binary" = "nb"),
        selected = c("m", "w", "nb")
      ),
      actionButton("go", "Analyze location"),
      helpText("Estimates are precomputed for all 50 states and DC, so results usually appear within a few seconds."),
      helpText(
        tags$a(href = "https://github.com/martensn/gaydar/tree/main", target = "_blank", "Git repository"),
        " · ",
        tags$a(href = "data.html", target = "_blank", "Data"),
        " · ",
        tags$a(href = "methods.html", target = "_blank", "Methodology")
      )
    ),

    mainPanel(
      div(
        class = "fill-height",

        div(
          class = "map-panel",
          shinycssloaders::withSpinner(leafletOutput("map", height = "100%"), color = "#cc8855")
        ),

        div(
          class = "donut-panel",
          shinycssloaders::withSpinner(plotlyOutput("composition_plot", height = "100%"), color = "#cc8855")
        )
      )
    )
  )
)

# ---------------------------
# Server
# ---------------------------

server <- function(input, output, session) {
  colors <- c(
    non_lgbt = "#57534a",
    lg       = "#cc8855",
    bi       = "#e3b23c",
    queer    = "#ff1493",
    trans    = "#b5567c"
  )
  label_map <- c(
    lgbt     = "LGBTQ",
    non_lgbt = "Straight",
    lg       = "Lesbian or Gay",
    bi       = "Bisexual",
    queer    = "Queer",
    trans    = "Transgender"
  )
  gender_labels <- c(m = "Men", w = "Women", nb = "Enbies")
  
  # ---- metric × gender → columns ----
  metric_gender_cols <- list(
    lgbt = list(m = "lgbt_m_map",  w = "lgbt_w_map",  nb = "lgbt_nb_map"),
    lg   = list(m = "lg_m_map",    w = "lg_w_map",    nb = "lg_nb_map"),
    bi   = list(m = "bi_m_map",    w = "bi_w_map",    nb = "bi_nb_map"),
    queer= list(m = "queer_m_map", w = "queer_w_map", nb = "queer_nb_map"),
    trans= list(m = "trans_m_map", w = "trans_w_map", nb = character(0))
  )
  
  get_identity_counts <- function(sf_obj, genders) {
    
    # Helper to sum gender-specific columns safely
    sum_cols <- function(cols) {
      if (length(cols) == 0) return(0)
      
      sf_obj |>
        sf::st_drop_geometry() |>
        dplyr::select(dplyr::all_of(cols)) |>
        dplyr::mutate(dplyr::across(everything(), as.numeric)) |>
        sum(na.rm = TRUE)
    }
    
    tibble::tibble(
      lg = sum_cols(unlist(metric_gender_cols$lg[genders])),
      bi = sum_cols(unlist(metric_gender_cols$bi[genders])),
      queer = sum_cols(unlist(metric_gender_cols$queer[genders])),
      trans = sum_cols(unlist(metric_gender_cols$trans[genders])),
      lgbt = sum_cols(unlist(metric_gender_cols$lgbt[genders]))
    )
  }
  
  # ---- geocode address ----
  address_point <- eventReactive(input$go, {
      
    req(input$street, input$city, input$state)
    
    addr_df <- tibble::tibble(
      street = input$street,
      city   = input$city,
      state  = input$state,
      zip    = ifelse(nzchar(input$zip), input$zip, NA_character_)
    )
    
    # Census geocoder
    geo <- tidygeocoder::geocode(
      addr_df,
      street  = "street",
      city    = "city",
      state   = "state",
      postalcode = "zip",
      method = "census",
      full_results = TRUE,
      api_options = list(census_return_type = "geographies")
    ) %>%
      mutate(
        county_fips = purrr::map_chr(
          geographies.Counties,
          ~ .x$GEOID[1]
        )
      ) %>%
      transmute(
        street,
        city,
        state,
        zip,
        lat,
        lon = long,
        county_fips
      )
    
    validate(
      need(!is.na(geo$lat) && !is.na(geo$lon),
           "Address could not be geocoded.")
    )
    
    geo
  })
  
  # ---- infer state ----
  selected_state <- reactive({
    geo <- address_point()
    toupper(geo$state)
  })
  
  point_sf <- reactive({
    geo <- address_point()
    
    st_as_sf(
      geo,
      coords = c("lon", "lat"),
      crs = 4326
    )
  })
  # ---- session-level state cache (avoids re-reading disk for repeated state) ----
  session_cache <- reactiveVal(list())

  # ---- load tract layer for inferred state ----
  bg_sf <- reactive({
    state <- selected_state()

    cache <- session_cache()
    if (!is.null(cache[[state]])) {
      message("Session cache hit: ", state)
      return(cache[[state]])
    }

    rates <- readRDS(file.path(root_dir,"data/hps/hps_acs_rates.rds"))

    cache_dir <- file.path(root_dir,"data/cache", "tract_state")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    cache_key  <- paste0("tract_", state, "_2023_calTRUE_g0.5.rds")
    cache_path <- file.path(cache_dir, cache_key)

    out <- if (file.exists(cache_path)) {
      message("Loading cached: ", cache_path)
      readRDS(cache_path)
    } else {
      message("Building state tract layer: ", state, " 2023")
      withProgress(
        message = paste0("Fetching live Census data for ", state, " (first request for a new state can take up to a minute)"),
        value = 0.2,
        {
          tmp <- build_tract_expected_layer(
            state_abbr      = state,
            year            = 2023,
            rates           = rates,
            use_calibration = TRUE,
            gamma           = 0.5
          )
          incProgress(0.5, detail = "Simplifying geometry")
          tmp <- tmp[sf::st_geometry_type(tmp) %in% c("POLYGON", "MULTIPOLYGON"), ]
          tmp <- rmapshaper::ms_simplify(tmp, keep = 0.05, keep_shapes = TRUE)
          tmp <- sf::st_make_valid(tmp)
          tmp <- tmp[!sf::st_is_empty(tmp), ]
          saveRDS(tmp, cache_path)
          tmp
        }
      )
    }

    # Ensure geometries are valid before caching — st_intersection downstream
    # will crash silently (caught by safe_summarize) on any bad geometry
    out <- sf::st_make_valid(out)
    out <- out[sf::st_is_valid(out) & !sf::st_is_empty(out), ]

    session_cache(modifyList(cache, setNames(list(out), state)))
    out
  })

  # ---- dissolve tracts to the selected map geography (tract/county/PUMA/CD) ----
  # Only the choropleth map itself respects this -- the radius/metro/state
  # donut comparisons below stay on the raw tract layer (bg_sf()), since
  # those depend on precise areal overlap (a 10-mile buffer, a CBSA) where
  # aggregating first to a coarser geography would badly overstate partial
  # overlaps (e.g. one huge rural PUMA barely clipped by the radius buffer
  # would otherwise contribute its entire population).
  geo_sf <- reactive({
    req(bg_sf(), input$geo_level)
    sf_obj <- bg_sf()
    level  <- input$geo_level

    if (level == "tract") return(sf_obj)

    if (level == "cd") {
      cd_lookup <- national_tract_rates |> dplyr::distinct(GEOID, cd_id)
      sf_obj <- sf_obj |> dplyr::left_join(cd_lookup, by = "GEOID")
    }

    dissolve_geo(sf_obj, geo_level_info[[level]]$col)
  })

  # ---- 10-mile radius ----
  radius_bg <- reactive({
    req(bg_sf(), point_sf(), input$radius_miles)
    
    pt <- st_transform(point_sf(), st_crs(bg_sf()))
    
    # convert miles → meters
    dist_m <- input$radius_miles * 1609.34
    
    buf <- st_buffer(pt, dist = dist_m)
    
    st_intersection(bg_sf(), buf)
  })
  metro_bg <- reactive({
    req(point_sf(), bg_sf())

    # Determine the user's CBSA by finding which tract in the (already
    # correctly-assigned, including spatial-fallback-resolved) tract layer
    # actually contains their point, rather than doing a second, independent
    # county-FIPS lookup against the static crosswalk CSV -- that separate
    # lookup used to be able to disagree with the tract layer's own cbsa_id
    # (e.g. for Connecticut, whose planning-region FIPS codes the crosswalk
    # doesn't map to a real CBSA), which was the source of the Connecticut
    # metro-area bug.
    pt <- sf::st_transform(point_sf(), sf::st_crs(bg_sf()))
    hit <- suppressWarnings(sf::st_join(pt, bg_sf() |> dplyr::select(cbsa_id), join = sf::st_intersects))
    cbsa_id <- hit$cbsa_id[1]

    validate(need(!is.na(cbsa_id), "Address not in a CBSA."))

    bg_sf() |> filter(cbsa_id == !!cbsa_id)
  })
  
  # ---- summarize area ----
  summarize_area <- function(sf_obj) {
    sf_obj |>
      st_drop_geometry() |>
      summarize(
        total_pop = sum(total_pop, na.rm = TRUE),
        
        lg    = sum(lg_m_map + lg_w_map + lg_nb_map, na.rm = TRUE),
        bi    = sum(bi_m_map + bi_w_map + bi_nb_map, na.rm = TRUE),
        queer = sum(queer_m_map + queer_w_map + queer_nb_map, na.rm = TRUE),
        trans = sum(trans_m_map + trans_w_map, na.rm = TRUE),
        lgbt = sum(lgbt_m_map + lgbt_w_map + lgbt_nb_map, na.rm=TRUE)
      ) |>
      mutate(non_lgbt= total_pop - lgbt)
  }
  
  # ---- combined summary ----
  area_summary <- reactive({
    bind_rows(
      summarize_area(radius_bg()) |>
        mutate(area = paste0(round(input$radius_miles,1), " Mile Radius")),
      summarize_area(metro_bg())  |> mutate(area = "Metro"),
      summarize_area(bg_sf())     |> mutate(area = "State")
    )
  })
  
  bg_map_layer <- reactive({
    req(geo_sf(), input$metric, input$gender, input$geo_level)

    sf_obj <- geo_sf()
    info   <- geo_level_info[[input$geo_level]]

    cols <- unlist(metric_gender_cols[[input$metric]][input$gender])
    req(length(cols) > 0)

    sf_obj <- sf_obj %>%
      mutate(
        numerator = rowSums(across(all_of(cols)), na.rm = TRUE),
        denominator = rowSums(
          across(all_of(paste0("total_", input$gender))),
          na.rm = TRUE
        ),
        shr_selected = numerator / pmax(denominator, 1)
      )

    # ---- state percentile: rank within this state's own geo_sf() rows ----
    has_pop <- sf_obj$denominator > 0
    state_ecdf <- if (any(has_pop)) stats::ecdf(sf_obj$shr_selected[has_pop]) else function(x) NA_real_
    sf_obj$state_pctile <- ifelse(has_pop, state_ecdf(sf_obj$shr_selected) * 100, NA_real_)

    # ---- national percentile: rank within the same metric/gender/geography
    # combo, computed nationally from national_tract_rates ----
    national_df <- aggregate_national(info$col) %>%
      mutate(
        nat_numerator   = rowSums(across(all_of(cols)), na.rm = TRUE),
        nat_denominator = rowSums(across(all_of(paste0("total_", input$gender))), na.rm = TRUE),
        nat_shr         = nat_numerator / pmax(nat_denominator, 1)
      ) %>%
      filter(nat_denominator > 0)

    national_ecdf <- if (nrow(national_df) > 0) stats::ecdf(national_df$nat_shr) else function(x) NA_real_
    sf_obj$national_pctile <- ifelse(has_pop, national_ecdf(sf_obj$shr_selected) * 100, NA_real_)

    # ---- human-readable name (county/PUMA/CD) instead of a raw FIPS-style
    # id; tracts have no such lookup and keep showing their GEOID ----
    if (!is.null(info$names)) {
      sf_obj <- sf_obj |> dplyr::left_join(geo_names[[info$names]], by = "GEOID")
      sf_obj$display_id <- dplyr::coalesce(sf_obj$name, sf_obj$GEOID)
    } else {
      sf_obj$display_id <- sf_obj$GEOID
    }

    ordinal_suffix <- function(n) {
      dplyr::case_when(
        n %% 100 %in% 11:13 ~ "th",
        n %% 10 == 1        ~ "st",
        n %% 10 == 2        ~ "nd",
        n %% 10 == 3        ~ "rd",
        TRUE                 ~ "th"
      )
    }
    # Value-cell HTML for a percentile row -- the word "percentile" is its
    # own span so it can be styled as an accent color separately from the number.
    fmt_pctile <- function(p) {
      r <- round(p)
      ifelse(
        is.na(p),
        "n/a (no population)",
        paste0(r, ordinal_suffix(r), " <span class=\"gd-pctile\">percentile</span>")
      )
    }

    tooltip_row <- function(label, value) {
      sprintf("<tr><td class=\"gd-label\">%s</td><td class=\"gd-value\">%s</td></tr>", label, value)
    }

    sf_obj %>%
      mutate(
        label = paste0(
          "<table class=\"gd-table\">",
          tooltip_row(info$label, display_id),
          tooltip_row("lgbtq+", formatC(numerator, format = "f", digits = 0, big.mark = ",")),
          tooltip_row("total", formatC(denominator, format = "f", digits = 0, big.mark = ",")),
          tooltip_row("percent", sprintf("%.1f%%", shr_selected * 100)),
          "</table>"
        ) |> lapply(htmltools::HTML),
        popup = paste0(
          "<table class=\"gd-table\">",
          tooltip_row(info$label, display_id),
          tooltip_row("lgbtq+", formatC(numerator, format = "f", digits = 0, big.mark = ",")),
          tooltip_row("total", formatC(denominator, format = "f", digits = 0, big.mark = ",")),
          tooltip_row("percent", sprintf("%.1f%%", shr_selected * 100)),
          tooltip_row("state", fmt_pctile(state_pctile)),
          tooltip_row("national", fmt_pctile(national_pctile)),
          "</table>"
        ) |> lapply(htmltools::HTML)
      )
  })
  
  # ---- composition plot ----
  output$composition_plot <- plotly::renderPlotly({
    
    # --- safe summaries (don’t let one failing area kill the whole plot) ---
    safe_summarize <- function(name, expr) {
      out <- tryCatch(expr, error = function(e) NULL)
      if (is.null(out)) return(NULL)
      out |> mutate(area = name)
    }
    
    # Build the three areas (radius label matches your area_summary())
    rad_name <- paste0(round(input$radius_miles, 1), " mile radius")

    s_rad   <- safe_summarize(rad_name, summarize_area(radius_bg()))
    s_metro <- safe_summarize("metro",  summarize_area(metro_bg()))
    s_state <- safe_summarize("state",  summarize_area(bg_sf()))

    sums <- dplyr::bind_rows(s_rad, s_metro, s_state)

    validate(need(nrow(sums) > 0, "Waiting for area summaries…"))

    df <- dplyr::bind_rows(lapply(seq_len(nrow(sums)), function(i) {
      row <- sums[i, ]
      area_name <- row$area
      # Pull the correct SF object for this area
      if (area_name == rad_name) {
        sf_obj <- radius_bg()
      } else if (area_name == "metro") {
        sf_obj <- metro_bg()
      } else if (area_name == "state") {
        sf_obj <- bg_sf()
      } else {
        stop("Unknown area: ", area_name)
      }
      # Identity counts for selected genders
      id_counts <- get_identity_counts(sf_obj, input$gender)
      # Denominator = selected-gender population
      denom <- sf_obj |>
        sf::st_drop_geometry() |>
        dplyr::select(dplyr::all_of(paste0("total_", input$gender))) |>
        dplyr::mutate(dplyr::across(everything(), as.numeric)) |>
        sum(na.rm = TRUE)
      
      id_counts |>
        mutate(
          non_lgbt = denom - lgbt,
          area = area_name
        ) |>
        select(area, non_lgbt, lg, bi, queer, trans)
      })) |>
      tidyr::pivot_longer(-area, names_to = "identity", values_to = "count") |>
      group_by(area) |>
      mutate(
        total = sum(count, na.rm = TRUE),
        prop  = ifelse(total > 0, count / total, NA_real_),
        hover = paste0(
          "<b>", label_map[identity], "</b><br>",
          "Count: ", scales::comma(count), "<br>",
          "Percent: ", scales::percent(prop, accuracy = 0.1))) |>
      ungroup()
    
    
    # Order areas the way you want (only those that exist)
    area_order <- c(rad_name, "metro", "state")
    areas <- area_order[area_order %in% unique(df$area)]
    validate(need(length(areas) >= 1, "Waiting for valid areas…"))
    
    pretty_gender <- paste(gender_labels[input$gender], collapse = ", ")
    
    # --- explicit domains for each donut ---
    n <- length(areas)
    x_starts <- seq(0, 1, length.out = n + 1)[1:n]
    x_ends   <- seq(0, 1, length.out = n + 1)[2:(n + 1)]
    x_mids   <- (x_starts + x_ends) / 2
    
    fig <- plotly::plot_ly()
    
    for (i in seq_along(areas)) {
      a <- areas[i]
      sub <- df |> filter(area == a)
      
      # One pie trace per area with its own x-domain
      fig <- fig |>
        plotly::add_trace(
          data = sub,
          type   = "pie",
          labels = ~label_map[identity],
          values = ~count,
          hole   = 0.6,
          textinfo  = "none",
          hoverinfo = "text",
          text = ~hover,
          marker = list(colors = unname(colors[as.character(sub$identity)])),
          domain = list(x = c(x_starts[i], x_ends[i]), y = c(0, 1)),
          showlegend = (i == 1)  # single legend
        )
    }
    
    # Center labels (paper coords so they land in each donut center)
    ann <- list()
    
    for (i in seq_along(areas)) {
      
      a <- areas[i]
      sub <- df |> filter(area == a)
      
      total_pop <- sum(sub$count)
      straight_pop <- sum(sub$count[sub$identity == "non_lgbt"], na.rm = TRUE)
      lgbt_pop <- total_pop - straight_pop
      
      # ---- CENTER TOTAL (inside donut) ----
      ann[[length(ann) + 1]] <- list(
        x = x_mids[i],
        y = 0.5,
        xref = "paper",
        yref = "paper",
        showarrow = FALSE,
        
        align   = "center",
        xanchor = "center",
        yanchor = "middle",
        text = paste0(
          "LGBTQ: ", scales::comma(lgbt_pop), "<br>",
          "Total: ", scales::comma(total_pop),
          "</span>"
        ),
        font = list(
          family = "Helvetica",
          size   = 14,
          color  = "#fffdda"
        )
      )

      # ---- GEOGRAPHY LABEL (above donut) ----
      ann[[length(ann) + 1]] <- list(
        x = x_mids[i],
        y = 1.05,
        xref = "paper",
        yref = "paper",
        showarrow = FALSE,
        
        align   = "center",
        xanchor = "center",
        yanchor = "bottom",
        
        text = paste0("<b>", a, "</b>"),
        
        font = list(
          family = "Helvetica",
          size   = 18,
          color  = "#cc8855"
        )
      )
    }

    fig |>
      plotly::layout(
        annotations = ann,
        margin = list(t = 70, b = 10, l = 10, r = 10),
        legend = list(orientation = "v", font = list(color = "#fffdda")),
        paper_bgcolor = "#282b33",
        plot_bgcolor  = "#282b33",

        font = list(
          family = "Helvetica",
          size   = 13,
          color  = "#fffdda"
        )
      )
  })

  output$map <- renderLeaflet({
    leaflet() |>
      addProviderTiles("CartoDB.DarkMatter")
  })
    
  # ---- map (simple highlight) ----
  observeEvent(input$go, {
    
    req(point_sf(), bg_map_layer())
    
    pt  <- st_transform(point_sf(), 4326)
    bgs <- st_transform(bg_map_layer(), 4326)

    # ---- FIX 1: geometry collections ----
    bgs <- bgs |>
      sf::st_make_valid()
    
    bgs <- bgs[!sf::st_is_empty(bgs), ]
    
    # keep only polygonal geometries
    geom_types <- as.character(sf::st_geometry_type(bgs))
    bgs <- bgs[geom_types %in% c("POLYGON", "MULTIPOLYGON"), ]
    
    # optional: harmonize type
    if (nrow(bgs) > 0) {
      bgs <- sf::st_cast(bgs, "MULTIPOLYGON", warn = FALSE)
    } 
    # ---- winsorized color scale ----
    # Tracts with zero selected-gender population have numerator = denominator
    # = 0, so shr_selected = 0 -- indistinguishable from a real, very-low
    # share at the dark end of the viridis scale. Excluded from the domain
    # calc and rendered as a flat neutral gray instead.
    vals <- bgs$shr_selected
    zero_pop <- !is.na(bgs$denominator) & bgs$denominator <= 0
    zero_pop_color <- "#d9d9d9"

    qs <- quantile(vals[!zero_pop], c(0.01, 0.99), na.rm = TRUE)
    vals_clamped <- pmin(pmax(vals, qs[1]), qs[2])

    # Warm gradient matching the site palette (muted dark -> orange -> pink)
    # instead of viridis's blue/green/yellow, which clashed with the theme.
    pal <- colorNumeric(
      palette = c("#3a3228", "#cc8855", "#ff1493"),
      domain = qs,
      na.color = "#00000000"
    )

    fill_colors <- ifelse(zero_pop, zero_pop_color, pal(vals_clamped))

    legend_title <- paste0(
      label_map[input$metric],
      " Percent<br>",
      "<span style='font-weight: normal;'>",
      paste(gender_labels[input$gender], collapse = ", "),
      "</span>"
    )
    
    leafletProxy("map") |>
      clearShapes() |>
      clearMarkers() |>
      clearControls() |>
      addPolygons(
        data = bgs,
        fillColor = fill_colors,
        fillOpacity = 0.7,
        color = "rgba(255, 253, 218, 0.35)",
        weight = 0.4,
        label = ~label,
        popup = ~popup,
        highlightOptions = highlightOptions(
          weight = 2,
          color = "#fffdda",
          bringToFront = TRUE
        )
      ) |>
      addCircles(
        lng = st_coordinates(pt)[1,1],
        lat = st_coordinates(pt)[1,2],
        radius = input$radius_miles * 1609.34,
        color = "#fffdda",
        weight = 1,
        fill = FALSE,
        opacity = 0.6
      ) |>
      addCircleMarkers(
        data = pt,
        radius = 6,
        color = "#ff1493",
        fillOpacity = 1
      ) |>
      addLegend(
        position = "bottomright",
        pal = pal,
        values = vals_clamped,
        title = legend_title,
        labFormat = labelFormat(
          transform = function(x) x * 100,
          suffix = "%",
          digits = 1
        )
      ) |>
      addLegend(
        position = "bottomright",
        colors = zero_pop_color,
        labels = "No population",
        opacity = 0.7
      ) |>
      setView(
        lng  = address_point()$lon,
        lat  = address_point()$lat,
        zoom = 11
      )
  })
  
  # ---- status ----
  #output$status <- renderPrint({
  #  list(
  #    address = paste(
  #      input$street,
  #      input$city,
  #      input$state,
  #      input$zip
  #    ),
  #    state   = selected_state(),
  #   summary = area_summary()
  #  )
  #})
}

shinyApp(ui, server)