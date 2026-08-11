# =============================================================================
# predict_fig_nonadjacent_v2.R
# -----------------------------------------------------------------------------
# Renders the PEHD illustration for NON-ADJACENT pairs of WWTPs.
# Uses the correct polygon-buffer semantics (same math as the paper),
# not the incorrect centroid-circle from the original figure.
#
# Two pairs are rendered:
#   - WWTP1 -> WWTP4   (30 km, moderate non-adjacent)
#   - WWTP1 -> WWTP15  (21 km, less extreme non-adjacent)
#
# Outputs (side-by-side, in the paper's parent folder):
#   figure_nonadjacent_v2.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(tidyverse); library(ggplot2); library(patchwork)
})

HERE <- "~/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW/Modelos AR con rho est. y rho=0"
setwd(HERE)
cat("Working from:", getwd(), "\n")

source(file.path(HERE, "R", "core.R"))
source(file.path(HERE, "R", "plot.R"))
source(file.path(HERE, "R", "analysis.R"))

# ============================================================================
# OVERRIDES: polygon-buffer semantics (same as v2 of the main figure)
# ============================================================================

captured_fraction <- function(polygon_geom, outer, dots, radius) {
  buf   <- st_sf(geometry = st_buffer(polygon_geom, dist = radius))
  inter <- st_intersection(buf, st_geometry(outer))
  hits  <- lengths(st_intersects(dots, inter)) > 0
  list(buffer = buf, intersection = inter, captured = hits,
       proportion = mean(hits))
}

find_capture_radius <- function(inner, outer, dots,
                                target_fraction = 0.5, ...) {
  inner_geom <- st_geometry(inner)
  dists_to_inner <- as.numeric(st_distance(dots, inner_geom))
  n <- length(dists_to_inner)
  target_k <- max(1L, as.integer(ceiling(target_fraction * n)))
  r_exact  <- sort.int(dists_to_inner, partial = target_k)[target_k]

  fin <- captured_fraction(inner_geom, outer, dots, r_exact)
  list(
    target_fraction = target_fraction,
    radius          = r_exact,
    proportion      = fin$proportion,
    buffer          = fin$buffer,
    intersection    = fin$intersection,
    captured        = fin$captured
  )
}

plot_capture <- function(inner, outer, dots, captured, buffer, intersection,
                         center_inner = NULL, center_outer = NULL,
                         title = "", subtitle = "",
                         colors = NULL, labels = NULL) {
  colors <- .merge_with_defaults(colors, default_colors())
  labels <- .merge_with_defaults(labels, default_labels())

  dots$captured <- captured

  ggplot() +
    geom_sf(data = outer, fill = colors$outer_fill,
            color = colors$outer_stroke, linewidth = 0.5) +
    geom_sf(data = inner, fill = colors$inner_fill,
            color = colors$inner_stroke, linewidth = 0.5) +
    geom_sf(data = intersection, fill = colors$intersection,
            color = NA, alpha = 0.45) +
    geom_sf(data = buffer, fill = NA, color = colors$buffer_stroke,
            linewidth = 0.5, linetype = "22") +
    geom_sf(data = dots, aes(color = captured), size = 0.22, alpha = 0.75) +
    scale_color_manual(
      values = c(`FALSE` = colors$dot_outside, `TRUE` = colors$dot_inside),
      labels = c(`FALSE` = "outside", `TRUE` = "inside"),
      name   = "Dot status"
    ) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}

# ============================================================================
# LOAD DATA
# ============================================================================
shp_batch       <- file.path(HERE, "data", "WWTP_Batch",
                             "WWTP_061621_mergedWWTPs to Batch_updatedon07192024.shp")
tracts_shp_path <- file.path(HERE, "data", "Census_Tracts_Population",
                             "WWTP Batch Service Areas_Census Tracts_Population.shp")

cat("Loading WWTP shapefile + tracts ...\n")
wwtp_shp_batch <- sf::read_sf(shp_batch) |>
  sf::st_transform(32615) |>
  mutate(BatchName = toupper(trimws(BatchName)))

batch_key_shp <- wwtp_shp_batch |>
  sf::st_drop_geometry() |>
  distinct(BatchName) |>
  arrange(BatchName) |>
  mutate(Batch = paste0("WWTP", row_number()))

wwtp_shp_batch <- wwtp_shp_batch |>
  left_join(batch_key_shp, by = "BatchName")

Batch_focus <- paste0("WWTP", 1:15)
wwtp_shp_focus <- wwtp_shp_batch |>
  filter(Batch %in% Batch_focus) |>
  arrange(match(Batch, Batch_focus))

tracts_raw <- sf::st_read(tracts_shp_path, quiet = TRUE) |>
  sf::st_make_valid() |>
  filter(sf::st_is_valid(geometry), !sf::st_is_empty(geometry))

tracts_valid <- tracts_raw |>
  distinct(GEOID, .keep_all = TRUE) |>
  sf::st_transform(32615)

cat(sprintf("  %d plants, %d unique tracts\n",
            nrow(wwtp_shp_focus), nrow(tracts_valid)))

# ============================================================================
# RENDER PAIR
# ============================================================================
render_pair_plot <- function(src_wwtp, tgt_wwtp) {
  inner <- wwtp_shp_focus |> filter(Batch == src_wwtp)
  outer <- wwtp_shp_focus |> filter(Batch == tgt_wwtp)

  res <- buffer_capture_analysis(
    inner           = inner,
    outer           = outer,
    tracts          = tracts_valid,
    proj_crs        = 32615,
    pop_col         = "POP24_5YR",
    labels          = list(inner = src_wwtp, outer = tgt_wwtp),
    target_fraction = 0.50,
    people_per_dot  = 200
  )

  ctarget <- find_capture_radius(inner, outer, res$dots,
                                 target_fraction = 0.50)

  plot_capture(
    inner, outer, res$dots,
    captured     = ctarget$captured,
    buffer       = ctarget$buffer,
    intersection = ctarget$intersection,
    title    = sprintf("%s → %s", src_wwtp, tgt_wwtp),
    subtitle = sprintf("r = %.2f km  |  achieved = %.3f",
                       ctarget$radius / 1000, ctarget$proportion)
  )
}

cat("\nRendering ADJACENT pair: WWTP12 -> WWTP13 (~1.6 km, share border) ...\n")
p_adj <- render_pair_plot("WWTP12", "WWTP13")

cat("Rendering NON-ADJACENT pair: WWTP1 -> WWTP3 (~17 km, no shared border) ...\n")
p_nonadj <- render_pair_plot("WWTP1", "WWTP3")

combined <- (p_adj | p_nonadj)

out_path <- file.path(HERE, "figure_pehd_adjacent_vs_nonadjacent_v2.pdf")
ggsave(out_path, combined, width = 12, height = 6)
cat("\n  -> ", out_path, "\n", sep = "")

cat("\nDone.\n")
