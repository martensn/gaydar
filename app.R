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
library(httr)
# jsonlite is used (via jsonlite::toJSON, fully-qualified below) but NOT
# attached with library() -- jsonlite::validate() masks shiny::validate(),
# which silently breaks every validate(need(...)) call in the app instead
# of erroring loudly, since jsonlite's version doesn't call stop().

options(tigris_use_cache = TRUE)

source("code/helpers.R")

# ---- travel-time isochrones (TravelTime API) ----
# Requires TRAVELTIME_ID / TRAVELTIME_KEY in .Renviron (free account at
# https://account.traveltime.com/, gitignored like CENSUS_API_KEY). Only
# used when the user picks "car" or "transit" for the radius mode --
# "distance" mode never touches the network, same as before.
#
# Returns a single-feature sf polygon (WGS84) reachable within `minutes` of
# (lat, lng) by `mode` ("driving" or "public_transport"), or throws an error
# with a message safe to surface via validate()/need().
fetch_isochrone <- function(lat, lng, minutes, mode) {
  app_id  <- Sys.getenv("TRAVELTIME_ID")
  api_key <- Sys.getenv("TRAVELTIME_KEY")
  if (!nzchar(app_id) || !nzchar(api_key)) {
    stop("Travel-time isochrones aren't configured on this deployment.")
  }

  body <- list(
    departure_searches = list(list(
      id = "iso",
      coords = list(lat = lat, lng = lng),
      departure_time = strftime(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      travel_time = minutes * 60,
      transportation = list(type = mode)
    ))
  )

  resp <- httr::POST(
    url = "https://api.traveltimeapp.com/v4/time-map",
    httr::add_headers(
      "Content-Type"     = "application/json",
      "Accept"           = "application/geo+json",
      "X-Application-Id" = app_id,
      "X-Api-Key"        = api_key
    ),
    body = jsonlite::toJSON(body, auto_unbox = TRUE)
  )

  if (httr::http_error(resp)) {
    if (httr::status_code(resp) == 429) {
      stop("Too many isochrone requests right now -- please try again in a minute.")
    }
    stop("Couldn't reach the isochrone service -- try again shortly.")
  }

  txt <- httr::content(resp, "text", encoding = "UTF-8")
  tmp <- tempfile(fileext = ".geojson")
  on.exit(unlink(tmp))
  writeLines(txt, tmp)
  poly <- sf::st_read(tmp, quiet = TRUE)

  if (nrow(poly) == 0) {
    stop("No reachable area found for that travel time.")
  }

  sf::st_make_valid(poly)
}

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

# ---- national county/PUMA/CD choropleth layers (built offline by
# code/build_national_geo_layers.R) -- the map renders these directly for
# those three geography levels instead of dissolving just the geocoded
# state, so panning covers the whole country. Tract stays per-state. ----
national_geo_layers <- readRDS(file.path(root_dir, "data/cache/national_geo_layers.rds"))

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
    tags$script(HTML("
      // Reconfigures the radius slider's range/step/unit when the mode
      // dropdown changes -- Shiny's updateSliderInput() can't touch the
      // ion.rangeSlider postfix, so this talks to the slider widget directly.
      var radiusModeConfig = {
        distance:          { min: 0, max: 50, step: 1, value: 10, postfix: ' mi' },
        driving:           { min: 5, max: 60, step: 5, value: 15, postfix: ' min' },
        public_transport:  { min: 5, max: 90, step: 5, value: 30, postfix: ' min' }
      };
      Shiny.addCustomMessageHandler('update_radius_slider', function(mode) {
        var cfg = radiusModeConfig[mode];
        var el = $('#radius_miles');
        var slider = el.data('ionRangeSlider');
        if (!cfg || !slider) return;
        slider.update({ min: cfg.min, max: cfg.max, step: cfg.step, from: cfg.value, postfix: cfg.postfix });
        el.trigger('change');
      });
    ")),
    tags$script(HTML("
      // Rainbow hover, matching the scintillate keyframe used for link
      // hovers on martensn.github.io -- applied here to the map's hovered
      // geography outline and the donut chart's hovered slice border,
      // since neither is a plain <a> text link CSS can animate directly.
      var GD_RAINBOW = ['#ff6b6b','#ffd93d','#6bcb77','#4d96ff','#9b59b6','#ff1493','#ff8c00'];

      // ---- map: delegated hover, attached to `document` rather than the
      // #map div -- this script runs in <head>, before Shiny has rendered
      // the leaflet widget into the body, so #map doesn't exist yet at
      // attach time. document always exists, and by the time a real hover
      // happens the tract paths certainly do too. Excludes the address-
      // point marker (also a path.leaflet-interactive) via its own class. ----
      (function() {
        var active = null, idx = 0, timer = null, orig = null;
        var SEL = 'path.leaflet-interactive:not(.gd-address-marker)';

        function start(path) {
          if (active === path) return;
          stop();
          active = path;
          orig = { stroke: path.getAttribute('stroke'), width: path.getAttribute('stroke-width') };
          path.setAttribute('stroke-width', 2.5);
          path.parentNode.appendChild(path); // bring to front
          idx = 0;
          timer = setInterval(function() {
            idx = (idx + 1) % GD_RAINBOW.length;
            path.setAttribute('stroke', GD_RAINBOW[idx]);
          }, 150);
        }
        function stop() {
          if (timer) clearInterval(timer);
          timer = null;
          if (active && orig) {
            active.setAttribute('stroke', orig.stroke);
            active.setAttribute('stroke-width', orig.width);
          }
          active = null;
        }
        document.addEventListener('mouseover', function(e) {
          var path = e.target.closest && e.target.closest(SEL);
          if (path) start(path);
        });
        document.addEventListener('mouseout', function(e) {
          var toEl = e.relatedTarget;
          var stillInside = toEl && toEl.closest && toEl.closest(SEL) === active;
          if (!stillInside) stop();
        });
      })();

      // ---- donut chart: cycle the hovered slice's border via Plotly's
      // JS API; re-attached (idempotently) on every Shiny render since
      // Plotly.react can occasionally hand back a fresh graph div ----
      function gdAttachDonutRainbow() {
        var gd = document.getElementById('composition_plot');
        if (!gd || !gd.on || gd._gdRainbowHooked) return;
        gd._gdRainbowHooked = true;
        var idx = 0, timer = null, trace = null, point = null, origColor = null;

        gd.on('plotly_hover', function(ev) {
          var pt = ev.points && ev.points[0];
          if (!pt || !pt.data || !pt.data.marker.line) return;
          trace = pt.curveNumber;
          point = pt.pointNumber;
          var lineColors = pt.data.marker.line.color.slice();
          origColor = lineColors[point];
          idx = 0;
          if (timer) clearInterval(timer);
          timer = setInterval(function() {
            idx = (idx + 1) % GD_RAINBOW.length;
            var colors = lineColors.slice();
            colors[point] = GD_RAINBOW[idx];
            Plotly.restyle(gd, { 'marker.line.color': [colors] }, [trace]);
          }, 150);
        });

        gd.on('plotly_unhover', function() {
          if (timer) clearInterval(timer);
          timer = null;
          if (trace !== null && point !== null) {
            var current = gd.data[trace].marker.line.color.slice();
            current[point] = origColor;
            Plotly.restyle(gd, { 'marker.line.color': [current] }, [trace]);
          }
          trace = null; point = null;
        });
      }
      $(document).on('shiny:value', function(e) {
        if (e.target && e.target.id === 'composition_plot') {
          setTimeout(gdAttachDonutRainbow, 0);
        }
      });
    ")),
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
  tags$title("Gaydar — LGBTQ Population Estimator"),
  div(
    class = "app-header",
    div(
      class = "app-title-group",
      tags$h2("Gaydar"),
      tags$nav(
        class = "app-links",
        tags$a(href = "https://github.com/martensn/gaydar/tree/main", target = "_blank", "repo"),
        " · ",
        tags$a(href = "data.html", target = "_blank", "data"),
        " · ",
        tags$a(href = "methods.html", target = "_blank", "methodology")
      )
    ),
    tags$a(
      href = "https://martensn.github.io",
      class = "home-arrow",
      `aria-label` = "back to homepage",
      "←"
    )
  ),
  sidebarLayout(
    sidebarPanel(
      textInput("street", "address", value = "3349 n halsted st", placeholder = "123 Main St"),
      textInput("city",   "city",    value = "chicago",           placeholder = "Chicago"),
      textInput("state",  "state",   value = "il",                placeholder = "IL"),
      textInput("zip",    "zip",     value = "60657",             placeholder = "60637"),
      selectInput(
        "radius_mode",
        "mode",
        choices = c(
          "distance" = "distance",
          "car"      = "driving",
          "transit"  = "public_transport"
        ),
        selected = "distance"
      ),
      sliderInput(
        "radius_miles",
        "radius",
        min   = 0,
        max   = 50,
        value = 10,
        step  = 1,
        ticks = FALSE,
        post  = " mi"
      ),
      selectInput(
        "geo_level",
        "geography",
        choices = c(
          "census tract"           = "tract",
          "county"                 = "county",
          "puma"                   = "puma",
          "congressional district" = "cd"
        ),
        selected = "tract"
      ),
      selectInput(
        "metric",
        "population",
        choices = c(
          "lgbt"            = "lgbt",
          "lesbian and gay" = "lg",
          "bisexual"        = "bi",
          "queer"           = "queer",
          "transgender"     = "trans"
        ),
        selected = "lgbt"
      ),

      checkboxGroupInput(
        "gender",
        "gender",
        choices = c("man" = "m", "woman" = "w", "non-binary" = "nb"),
        selected = c("m", "w", "nb")
      ),
      actionButton("go", "analyze")
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
          tags$button(
            class = "donut-toggle",
            `aria-label` = "toggle charts",
            onclick = "
              var p = this.closest('.donut-panel');
              p.classList.toggle('collapsed');
              this.textContent = p.classList.contains('collapsed') ? '▲' : '▼';
              setTimeout(function() {
                var widget = HTMLWidgets.find('#map');
                if (widget && widget.getMap) widget.getMap().invalidateSize();
              }, 220);
            ",
            "▼"
          ),
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

  # ---- reconfigure the radius slider's range/unit when mode changes ----
  observeEvent(input$radius_mode, {
    session$sendCustomMessage("update_radius_slider", input$radius_mode)
  }, ignoreInit = TRUE)

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

  # ---- map geography (tract/county/PUMA/CD) ----
  # Only the choropleth map itself respects this -- the radius/metro/state
  # donut comparisons below stay on the raw tract layer (bg_sf()), since
  # those depend on precise areal overlap (a 10-mile buffer, a CBSA) where
  # aggregating first to a coarser geography would badly overstate partial
  # overlaps (e.g. one huge rural PUMA barely clipped by the radius buffer
  # would otherwise contribute its entire population).
  #
  # Tract stays scoped to the geocoded state (bg_sf()); county/PUMA/CD use
  # the offline-precomputed *national* layers (code/build_national_geo_layers.R)
  # instead of dissolving just the current state, so the choropleth covers
  # the whole country and users can pan beyond the searched address without
  # re-geocoding.
  geo_sf <- reactive({
    req(input$geo_level)
    level <- input$geo_level

    if (level == "tract") {
      req(bg_sf())
      return(bg_sf())
    }

    national_geo_layers[[level]]
  })

  # ---- travel-time isochrone: only (re-)fetched on "analyze", never on a
  # bare slider drag -- the API is rate-limited, unlike the free circular
  # buffer below, which can recompute freely as sliders move. ----
  isochrone_shape <- eventReactive(input$go, {
    req(input$radius_mode != "distance")

    geo <- address_point()
    iso <- tryCatch(
      fetch_isochrone(geo$lat, geo$lon, input$radius_miles, input$radius_mode),
      error = function(e) e
    )
    if (inherits(iso, "error")) {
      validate(need(FALSE, conditionMessage(iso)))
    }
    iso
  }, ignoreInit = TRUE)

  # ---- radius shape: circular buffer (distance mode) or the travel-time
  # isochrone fetched above (car/transit mode), in bg_sf()'s CRS ----
  radius_shape <- reactive({
    req(bg_sf(), point_sf(), input$radius_miles, input$radius_mode)

    pt <- st_transform(point_sf(), st_crs(bg_sf()))

    if (input$radius_mode == "distance") {
      dist_m <- input$radius_miles * 1609.34
      return(st_buffer(pt, dist = dist_m))
    }

    st_transform(isochrone_shape(), st_crs(bg_sf()))
  })

  # ---- area within the radius/isochrone ----
  radius_bg <- reactive({
    req(bg_sf(), radius_shape())
    st_intersection(bg_sf(), radius_shape())
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

    # ---- state percentile ----
    # Tract level: geo_sf() is already scoped to the geocoded state, so
    # rank directly within it. County/PUMA/CD: geo_sf() is now the full
    # national layer (so the map can be panned anywhere), so "state" has
    # to mean whichever state *that specific row* is in -- not the
    # originally-searched address's state -- via a grouped ecdf keyed by
    # each GEOID's own state-FIPS prefix.
    has_pop <- sf_obj$denominator > 0

    if (input$geo_level == "tract") {
      state_ecdf <- if (any(has_pop)) stats::ecdf(sf_obj$shr_selected[has_pop]) else function(x) NA_real_
      sf_obj$state_pctile <- ifelse(has_pop, state_ecdf(sf_obj$shr_selected) * 100, NA_real_)
    } else {
      sf_obj$row_state_fips <- substr(sf_obj$GEOID, 1, 2)
      pctile_lookup <- sf_obj %>%
        sf::st_drop_geometry() %>%
        filter(denominator > 0) %>%
        group_by(row_state_fips) %>%
        mutate(state_pctile = { e <- stats::ecdf(shr_selected); e(shr_selected) * 100 }) %>%
        ungroup() %>%
        select(GEOID, state_pctile)
      sf_obj <- sf_obj %>%
        left_join(pctile_lookup, by = "GEOID") %>%
        select(-row_state_fips)
    }

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
          tooltip_row("share", sprintf("%.1f percent", shr_selected * 100)),
          "</table>"
        ) |> lapply(htmltools::HTML),
        popup = paste0(
          "<table class=\"gd-table\">",
          tooltip_row(info$label, display_id),
          tooltip_row("lgbtq+", formatC(numerator, format = "f", digits = 0, big.mark = ",")),
          tooltip_row("total", formatC(denominator, format = "f", digits = 0, big.mark = ",")),
          tooltip_row("share", sprintf("%.1f percent", shr_selected * 100)),
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
    rad_name <- switch(input$radius_mode,
      distance         = paste0(round(input$radius_miles, 1), " mile radius"),
      driving          = paste0(round(input$radius_miles, 1), " min drive"),
      public_transport = paste0(round(input$radius_miles, 1), " min transit")
    )

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
          marker = list(
            colors = unname(colors[as.character(sub$identity)]),
            # baseline border matches the panel background (invisible at
            # rest); JS cycles a hovered slice's line color through the
            # site's rainbow sequence, same as link/button hovers
            line = list(color = rep("#32353f", nrow(sub)), width = 2)
          ),
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

    req(point_sf(), bg_map_layer(), radius_shape())

    pt  <- st_transform(point_sf(), 4326)
    bgs <- st_transform(bg_map_layer(), 4326)
    radius_outline <- st_transform(radius_shape(), 4326)

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

    # Warm gradient matching the site palette (muted khaki -> orange -> pink)
    # instead of viridis's blue/green/yellow, which clashed with the theme.
    # Low end is a warm khaki leaning toward the site's cream text color,
    # without going so light the map reads as just shades of pink.
    map_gradient <- c("#c5c0a5", "#cc8855", "#ff1493")
    pal <- colorNumeric(
      palette = map_gradient,
      domain = qs,
      na.color = "#00000000"
    )

    fill_colors <- ifelse(zero_pop, zero_pop_color, pal(vals_clamped))

    legend_title <- paste0(
      tolower(label_map[input$metric]), " percent<br>",
      "<span style='font-weight: normal;'>",
      tolower(paste(gender_labels[input$gender], collapse = ", ")),
      "</span>"
    )

    # ---- single-box legend: gradient bar + "no population" swatch in one
    # control, rather than two separately-positioned addLegend() boxes ----
    legend_ticks <- seq(qs[1], qs[2], length.out = 5)
    tick_positions <- seq(0, 100, length.out = 5)
    tick_html <- paste0(
      mapply(function(v, pos) {
        sprintf('<span class="gd-legend-tick" style="bottom:%.2f%%">%.1f%%</span>', pos, v * 100)
      }, legend_ticks, tick_positions),
      collapse = ""
    )
    legend_html <- sprintf(
      '<div class="gd-legend-title">%s</div>
       <div class="gd-legend-gradient" style="background: linear-gradient(to top, %s, %s, %s);">%s</div>
       <div class="gd-legend-nopop"><span class="gd-legend-swatch"></span>no population</div>',
      legend_title, map_gradient[1], map_gradient[2], map_gradient[3], tick_html
    )

    leafletProxy("map") |>
      clearShapes() |>
      clearMarkers() |>
      clearControls() |>
      addPolygons(
        data = bgs,
        # layerId + a stable class so the client-side rainbow-hover script
        # (see tags$script above) can target these paths via CSS/DOM only
        # -- highlightOptions is deliberately omitted; JS owns the hover
        # styling here so it can animate through the rainbow sequence
        # instead of a single static highlight color.
        fillColor = fill_colors,
        fillOpacity = 0.7,
        color = "rgba(255, 253, 218, 0.35)",
        weight = 0.4,
        label = ~label,
        popup = ~popup
      ) |>
      addPolygons(
        data = radius_outline,
        options = pathOptions(interactive = FALSE),
        color = "#fffdda",
        weight = 1,
        fill = FALSE,
        opacity = 0.6
      ) |>
      addCircleMarkers(
        data = pt,
        radius = 6,
        color = "#ff1493",
        fillOpacity = 1,
        options = pathOptions(className = "gd-address-marker")
      ) |>
      addControl(
        html = HTML(legend_html),
        position = "bottomright"
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