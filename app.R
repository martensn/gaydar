# app.R -- primitive Shiny prototype for BG-level LGB mapping
# ---------------------------------------------------------------
# Assumptions:
#   - You have precomputed/serialized HPS-derived rates as an RDS at:
#       data/hps_acs_rates.rds
#     with columns: state_abbr, sex ("M"/"F"), acs_bin, lg_shr, lgb_shr
#   - You have a Census API key configured for tidycensus:
#       tidycensus::census_api_key("YOUR_KEY", install = TRUE)
#
# Minimal workflow:
#   - Choose a state and ACS year
#   - App pulls ACS B01001 BG counts + geometry via tidycensus
#   - Applies demographic rates to generate expected LG/LGB per BG
#   - Optional: county-level calibration using same-sex couples proxy
#
# Notes:
#   - This is intentionally a FIRST PASS to get mapping working.
#   - BG mapping is heavy; we cache per-state/year results to disk.

library(shiny)
library(leaflet)
library(sf)
library(dplyr)
library(tidyr)
library(stringr)
library(tidycensus)
library(tigris)

options(tigris_use_cache = TRUE)

source("code/helpers.R")

ui <- fluidPage(
  titlePanel("Gaydar"),
  sidebarLayout(
    sidebarPanel(
      selectInput("state_abbr", "State", choices = get_state_choices(), selected = "IL"),
      numericInput("year", "ACS Vintage", value = 2023, min = 2010, max = 2024, step = 1),
      #checkboxInput("use_calibration", "Calibrate based on geography of same-sex couples?", value = TRUE),
      #sliderInput("gamma", "Calibration strength (gamma)", min = 0, max = 2, value = 0.5, step = 0.1),
      selectInput(
        "metric",
        "Population:",
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
        "Gender",
        choices = c("Man"="m","Woman"="w","Non-binary"="nb"),
        selected = c("m","w","nb"),
        inline = FALSE
      )
      #checkboxInput("simplify", "Simplify geometry (faster)", value = TRUE),
      #sliderInput("simplify_tol", "Simplify tolerance (degrees)", min = 0.0001, max = 0.01, value = 0.002, step = 0.0001)
    ),
    mainPanel(
      leafletOutput("map", height = 700),
      tags$hr(),
      verbatimTextOutput("status")
    )
  )
)

server <- function(input, output, session) {
  metric_gender_cols <- list(
    lgbt = list(
      m  = "lgbt_m_map",
      w  = "lgbt_w_map",
      nb = "lgbt_nb_map"
    ),
    lg = list(
      m  = "lg_m_map",
      w  = "lg_w_map",
      nb = "nb_nt_lg"
    ),
    bi = list(
      m  = "bi_m_map",
      w  = "bi_w_map",
      nb = "bi_nb_map"
    ),
    queer = list(
      m  = "queer_m_map",
      w  = "queer_w_map",
      nb = "queer_nb_map"
    ),
    trans = list(
      m  = "trans_m_map",
      w  = "trans_w_map",
      nb = character(0)
    )
  )
  
  # base leaflet
  output$map <- renderLeaflet({
    leaflet() |>
      addProviderTiles("CartoDB.Positron") |>
      setView(lng = -89.5, lat = 39.8, zoom = 6)
  })

  bg_sf <- reactive({
    validate(
      need(file.exists("data/hps/hps_acs_rates.rds"),
           "Missing data/hps_acs_rates.rds. Create it from your HPS pipeline first.")
    )

    rates <- readRDS("data/hps/hps_acs_rates.rds")
    state <- input$state_abbr
    year  <- input$year

    # cache on disk to avoid repeated heavy ACS pulls
    cache_dir <- file.path("data/cache", "bg_state")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    cache_key <- paste0("bg_", state, "_", year, "_calTRUE_g0.5.rds")
    cache_path <- file.path(cache_dir, cache_key)

    if (file.exists(cache_path)) {
      message("Loading cached: ", cache_path)
      return(readRDS(cache_path))
    }

    message("Building state BG layer: ", state, " ", year)

    # build expected counts/shares at BG level
    out <- build_bg_expected_layer(
      state_abbr = state,
      year = year,
      rates = rates,
      use_calibration = input$use_calibration,
      gamma = 0.5
    )
    
    saveRDS(out, cache_path)
    out
  })
  sf_view <- reactive({
    
    sf_obj <- bg_sf()
    req(sf_obj, input$metric, input$gender)
    
    validate(
      need(length(input$gender) > 0, "Select at least one gender.")
    )
    
    # ---- map (metric × gender) → columns ----
    cols <- unlist(metric_gender_cols[[input$metric]][input$gender])
    
    validate(
      need(length(cols) > 0, "No columns for this selection.")
    )
    
    sf_obj %>%
      mutate(
        est_selected = rowSums(across(all_of(cols)), na.rm = TRUE),
        denom_selected = rowSums(across(all_of(paste0("total_", input$gender))), na.rm = TRUE),
        shr_selected = est_selected / pmax(denom_selected, 1)
      )
  })
    
  observe({
    
    sf_obj <- sf_view()
    cat("metric:", input$metric,
        "gender:", paste(input$gender, collapse=","),
        "range:", paste(range(sf_obj$shr_selected, na.rm=TRUE), collapse=" to "),
        "\n")
    req(sf_obj)
    
    # ---- geometry prep (unchanged) ----
    sf_obj <- sf_obj |> sf::st_make_valid()
    
    geom_type <- sf::st_geometry_type(sf_obj, by_geometry = TRUE)
    sf_obj <- sf_obj[geom_type %in% c("POLYGON", "MULTIPOLYGON"), ]
    
    sf_obj <- suppressWarnings(
      sf::st_simplify(
        sf_obj,
        dTolerance = 0.002,
        preserveTopology = TRUE
      )
    )
    
    
    # Ensure leaflet-safe geometry
    sf_obj <- sf_obj |>
      sf::st_cast("MULTIPOLYGON", warn = FALSE)
    
    # Stable domain based on ALL genders for the current metric
    vals <- sf_obj$shr_selected
    qs <- quantile(vals, probs = c(0.005, 0.995), na.rm = TRUE)
    vals_clamped <- pmin(pmax(vals, qs[1]), qs[2])
    
    pal <- colorNumeric("viridis", domain = qs, na.color = "#00000000")

    lbl <- sprintf(
      "<strong>GEOID:</strong> %s<br/>
   Pop (18+): %s<br/>
   E[selected]: %s<br/>
   Share: %s",
      sf_obj$GEOID,
      format(round(sf_obj$total_pop), big.mark = ","),
      format(round(sf_obj$est_selected), big.mark = ","),
      paste0(formatC(100 * sf_obj$shr_selected, digits = 0, format = "f"), "%")
    ) |> lapply(htmltools::HTML)

    leafletProxy("map", data = sf_obj) |>
      clearShapes() |>
      clearControls() |>
      addPolygons(
        weight = 0.2,
        opacity = 1,
        color = "#666666",
        fillOpacity = 0.7,
        fillColor = ~pal(vals_clamped),
        label = lbl,
        highlightOptions = highlightOptions(weight = 1.2, color = "#444444", bringToFront = TRUE)
      ) |>
      addLegend(
        position = "bottomright",
        pal = pal,
        values = vals_clamped,
        title = paste0(input$metric," share among ",paste(input$gender, collapse = ", ")),
        opacity = 0.7
      )

    # zoom to state bounding box
    bb <- sf::st_bbox(sf_obj)
    leafletProxy("map") |>
      fitBounds(bb["xmin"], bb["ymin"], bb["xmax"], bb["ymax"])
  })

  output$status <- renderPrint({
    sf_obj <- sf_view()
    req(sf_obj)
    
    tibble::tibble(
      state = input$state_abbr,
      year = input$year,
      metric = input$metric,
      genders = paste(input$gender, collapse = ","),
      n_bg = nrow(sf_obj),
      total_pop = sum(sf_obj$total_pop, na.rm = TRUE),
      est_selected = sum(sf_obj$est_selected, na.rm = TRUE),
      shr_selected = sum(sf_obj$est_selected, na.rm = TRUE) / sum(sf_obj$total_pop, na.rm = TRUE)
    )
  })
}

shinyApp(ui, server)
