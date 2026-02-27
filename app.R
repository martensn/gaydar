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
      numericInput("year", "ACS year (5-year)", value = 2023, min = 2010, max = 2024, step = 1),
      checkboxInput("use_calibration", "Calibrate based on geography of same-sex couples?", value = TRUE),
      sliderInput("gamma", "Calibration strength (gamma)", min = 0, max = 2, value = 0.5, step = 0.1),
      selectInput(
        "metric",
        "Map metric",
        choices = c(
          "LGBT Share" = "shr_lgbt_map",
          "Lesbian and Gay Share" = "shr_lg_map",
          "Bisexual Share" = "shr_bi_map",
          "Queer Share" = "shr_queer_map",
          "Trans Share" = "shr_trans_map"
        ),
        selected = "shr_lgbt_map"
      ),
      checkboxInput("simplify", "Simplify geometry (faster)", value = TRUE),
      sliderInput("simplify_tol", "Simplify tolerance (degrees)", min = 0.0001, max = 0.01, value = 0.002, step = 0.0001)
    ),
    mainPanel(
      leafletOutput("map", height = 700),
      tags$hr(),
      verbatimTextOutput("status")
    )
  )
)

server <- function(input, output, session) {

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

    cache_key <- paste0("bg_", state, "_", year, "_cal", input$use_calibration, "_g", input$gamma, ".rds")
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
      gamma = input$gamma
    )
    
    saveRDS(out, cache_path)
    out
  })

  observeEvent(bg_sf(), {
    sf_obj <- bg_sf()

    # ---- ALWAYS normalize geometry for leaflet ----
    sf_obj <- sf_obj |>
      sf::st_make_valid()
    
    # Drop geometry collections safely (keep polygons + multipolygons)
    geom_type <- sf::st_geometry_type(sf_obj, by_geometry = TRUE)
    
    sf_obj <- sf_obj[geom_type %in% c("POLYGON", "MULTIPOLYGON"), ]
    
    # optional simplification for browser speed
    if (isTRUE(input$simplify)) {
      sf_obj <- suppressWarnings(
        sf::st_simplify(
          sf_obj,
          dTolerance = input$simplify_tol,
          preserveTopology = TRUE
        )
      )
    }
    
    # Ensure leaflet-safe geometry
    sf_obj <- sf_obj |>
      sf::st_cast("MULTIPOLYGON", warn = FALSE)
    
    # Final safety check
    if (nrow(sf_obj) == 0) {
      stop("All geometries dropped after normalization — unexpected geometry types.")
    }

    metric <- input$metric
    vals <- sf_obj[[metric]]
    dom <- vals[is.finite(vals)]
    
    validate(
      need(length(dom) > 0, paste0("Metric '", metric, "' has no finite values to plot."))
    )
    
    pal <- colorNumeric("viridis", domain = dom, na.color = "#00000000")

    lbl <- sprintf(
      "<strong>GEOID:</strong> %s<br/>Pop (18+): %s<br/>E[LGBT] (uncal): %s<br/>E[LGBT] (calib): %s<br/>Share: %s",
      sf_obj$GEOID,
      format(round(sf_obj$pop18), big.mark = ",", scientific=FALSE),
      format(round(sf_obj$est_lgbt_uncal), big.mark = ",", scientific=FALSE),
      format(round(sf_obj$est_lgbt_map), big.mark = ",", scientific=FALSE),
      paste0(formatC(100 * sf_obj$shr_lgbt_map, format = "f", digits = 1), "%")
    ) |> lapply(htmltools::HTML)

    leafletProxy("map", data = sf_obj) |>
      clearShapes() |>
      clearControls() |>
      addPolygons(
        weight = 0.2,
        opacity = 1,
        color = "#666666",
        fillOpacity = 0.7,
        fillColor = ~pal(vals),
        label = lbl,
        highlightOptions = highlightOptions(weight = 1.2, color = "#444444", bringToFront = TRUE)
      ) |>
      addLegend(
        position = "bottomright",
        pal = pal,
        values = vals,
        title = metric,
        opacity = 0.7
      )

    # zoom to state bounding box
    bb <- sf::st_bbox(sf_obj)
    leafletProxy("map") |>
      fitBounds(bb["xmin"], bb["ymin"], bb["xmax"], bb["ymax"])
  })

  output$status <- renderPrint({
    req(bg_sf())
    sf_obj <- bg_sf()
    
    # ---- state-level status ----
    state_status <- tibble::tibble(
      state = input$state_abbr,
      year = input$year,
      n_bg = nrow(sf_obj),
      pop18_total = sum(sf_obj$pop18, na.rm = TRUE),
      est_lgbt_uncal = sum(sf_obj$est_lgbt_uncal, na.rm = TRUE),
      est_lgbt_map   = sum(sf_obj$est_lgbt_map, na.rm = TRUE),
      share_lgbt_map =
        sum(sf_obj$est_lgbt_map, na.rm = TRUE) /
        sum(sf_obj$pop18, na.rm = TRUE)
    )
    
    # ---- metro-level summary ----
    metro_status <- sf_obj %>%
      mutate(
        cbsa_name = if_else(
          is.na(cbsa_name),
          paste(input$state_abbr, "non-CBSA"),
          cbsa_name
        )
      ) %>%
      st_drop_geometry() %>%
      group_by(cbsa_name) %>%
      summarise(
        pop18 = sum(pop18, na.rm = TRUE),
        est_lgbt_uncal = sum(est_lgbt_uncal, na.rm = TRUE),
        est_lgbt_map   = sum(est_lgbt_map, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        shr_lgbt_uncal = est_lgbt_uncal / pop18,
        shr_lgbt_map   = est_lgbt_map   / pop18
      ) %>%
      arrange(desc(pop18))
    
    # ---- printed output ----
    list(
      state_summary = state_status,
      metro_summary = metro_status
    )
  })
}

shinyApp(ui, server)
