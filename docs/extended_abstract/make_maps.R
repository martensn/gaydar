# docs/extended_abstract/make_maps.R
# One composite page: three metro tract choropleths side by side, shared colour
# scale, shared legend. Boundaries are the full-resolution cartographic files
# from tigris (1:500k, clipped to shoreline), not the 5%-simplified geometry in
# the caches. Estimates are joined onto them by GEOID.
root_dir <- "."
suppressPackageStartupMessages({library(sf); library(dplyr); library(purrr)
                                library(ggplot2); library(scales); library(cowplot); library(tigris)})
options(tigris_use_cache = TRUE); sf::sf_use_s2(FALSE)
CACHE <- "data/cache/tract_state_agenid"

# Tighter, taller boxes on each urban core so three fit across a landscape page.
# All three boxes at the same aspect (height/width = 1.35 in projected terms)
# so the panels come out the same height and their titles align.
cities <- list(
  list(name="New York",    states=c("NY","NJ"), bb=c(-74.10, 40.56, -73.75, 40.92)),
  list(name="Los Angeles", states="CA",         bb=c(-118.60, 33.73, -118.10, 34.29)),
  list(name="Chicago",     states="IL",         bb=c(-87.95, 41.63, -87.52, 42.06))
)

# estimates at tract grain (sum the PUMA pieces), geometry dropped
est <- function(states) map_dfr(states, function(st)
    sf::st_drop_geometry(readRDS(file.path(CACHE, sprintf("tract_%s_2023_calTRUE_g0.5.rds", st))))) |>
  mutate(lgbt = lgbt_m_map + lgbt_w_map + lgbt_nb_map) |>
  group_by(GEOID) |> summarise(pop=sum(total_pop, na.rm=TRUE), lgbt=sum(lgbt, na.rm=TRUE), .groups="drop") |>
  filter(pop > 0) |> mutate(share = 100*lgbt/pop)

# full-resolution boundaries
geo <- function(states) map_dfr(states, ~ tigris::tracts(state=.x, year=2023, cb=TRUE)) |>
  select(GEOID) |> sf::st_transform(4326)

load_city <- function(ct) {
  bb <- sf::st_bbox(c(xmin=ct$bb[1], ymin=ct$bb[2], xmax=ct$bb[3], ymax=ct$bb[4]), crs=4326)
  g  <- suppressWarnings(sf::st_crop(geo(ct$states), bb))
  inner_join(g, est(ct$states), by="GEOID") |> mutate(city = ct$name)
}
dat <- map(cities, load_city)
allshare <- unlist(map(dat, ~ .x$share))
lims <- c(quantile(allshare,.02), quantile(allshare,.98))
cat(sprintf("shared scale %.1f%%-%.1f%%\n", lims[1], lims[2]))
walk(dat, ~ cat(sprintf("  %-12s %4d tracts  %.1f%%-%.1f%%  median %.1f%%\n",
  .x$city[1], nrow(.x), min(.x$share), max(.x$share), median(.x$share))))

pal <- c("#c5c0a5", "#cc8855", "#ff1493")
# patchwork 1.3.0's guide collection needs ggplot2 >= 3.5; this machine has
# 3.4.2. So: draw the panels legend-free, pull one legend out with cowplot, and
# lay it under all three as its own row.
panel <- function(d, bb, legend = FALSE) ggplot(d) +
  geom_sf(aes(fill = share), colour = "white", linewidth = 0.06) +
  scale_fill_gradientn(colours = pal, limits = lims, oob = scales::squish,
                       labels = function(x) paste0(x, "%"), name = "LGBTQ share of adults") +
  # lock the panel to its bbox rather than the data extent, so dropped water
  # tracts at the edges cannot change the panel's height
  coord_sf(xlim = c(bb[1], bb[3]), ylim = c(bb[2], bb[4]),
           expand = FALSE, datum = NA) +
  labs(title = d$city[1]) +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face="bold", size=17, hjust=0, margin=margin(b=4)),
        plot.margin = margin(2, 6, 2, 6),
        legend.position = if (legend) "bottom" else "none",
        legend.key.width = unit(2.8, "cm"), legend.key.height = unit(0.45, "cm"),
        legend.title = element_text(size=11, vjust=0.85),
        plot.background = element_rect(fill="white", colour=NA))

BB <- lapply(cities, `[[`, "bb")
leg  <- cowplot::get_legend(panel(dat[[1]], BB[[1]], legend = TRUE))
maps <- cowplot::plot_grid(panel(dat[[1]], BB[[1]]), panel(dat[[2]], BB[[2]]), panel(dat[[3]], BB[[3]]),
                           nrow = 1, align = "h", axis = "tb")
p <- cowplot::plot_grid(maps, leg, ncol = 1, rel_heights = c(1, 0.07)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("docs/extended_abstract/figs/maps_three_cities.pdf", p,
       width = 10.2, height = 7.0, device = cairo_pdf, bg = "white")
cat("wrote docs/extended_abstract/figs/maps_three_cities.pdf\n")
