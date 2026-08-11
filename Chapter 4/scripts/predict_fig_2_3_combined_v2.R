# =============================================================================
# predict_fig_2_3_combined_v2.R
# -----------------------------------------------------------------------------
# Regenerates the PEHD illustration using a buffer around the POLYGON boundary
# (mathematically correct) rather than a circle around the CENTROID
# (previous / incorrect illustration).
#
# The math in the paper says:
#     dist(x, Omega_i) = inf_{y in Omega_i} ||x - y||_2
# So the "buffer of radius r" grows from the BOUNDARY of Omega_i outward,
# not from a single centroid point.
#
# Output:
#   figure_2_3_combined_v2.pdf   (polygon-buffer version)
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
# OVERRIDES: use polygon-buffer semantics instead of centroid-circle
# ============================================================================

# Fraction of dots that fall inside (buffer(polygon, r) ∩ outer)
captured_fraction <- function(polygon_geom, outer, dots, radius) {
  buf   <- st_sf(geometry = st_buffer(polygon_geom, dist = radius))
  inter <- st_intersection(buf, st_geometry(outer))
  hits  <- lengths(st_intersects(dots, inter)) > 0
  list(
    buffer       = buf,
    intersection = inter,
    captured     = hits,
    proportion   = mean(hits)
  )
}

# Find smallest r such that buffer(Omega_i, r) captures target_fraction
# of the population dots inside Omega_j. Uses distances from dots to the
# INNER POLYGON (not to its centroid).
find_capture_radius <- function(inner, outer, dots,
                                target_fraction = 0.5,
                                frac_tol        = 0.005,
                                radius_tol_m    = 25,
                                max_iters       = 60) {
  inner_geom <- st_geometry(inner)

  # Precompute distances from every dot to the inner polygon (0 if inside)
  dists_to_inner <- as.numeric(st_distance(dots, inner_geom))
  r_hi <- max(dists_to_inner, na.rm = TRUE) * 1.05
  r_lo <- 0

  # Fast fraction check: sort distances, find k-th smallest s.t. k/n >= target
  n <- length(dists_to_inner)
  target_k <- max(1L, as.integer(ceiling(target_fraction * n)))
  r_exact  <- sort.int(dists_to_inner, partial = target_k)[target_k]

  # Sanity: bisect for consistency with the previous API (same-signature output)
  cap_at <- function(r) captured_fraction(inner_geom, outer, dots, r)

  # Use r_exact as our answer; the buffer/intersection are computed once at r_exact
  r_mid    <- r_exact
  fin      <- cap_at(r_mid)
  frac_mid <- fin$proportion

  list(
    target_fraction = target_fraction,
    radius          = r_mid,
    proportion      = frac_mid,
    iterations      = 1L,
    iter_log        = data.frame(iter = 1, r_lo = 0, r_hi = r_hi,
                                 r_mid = r_mid, frac_mid = frac_mid),
    buffer          = fin$buffer,
    intersection    = fin$intersection,
    captured        = fin$captured
  )
}

# Redraw plot_capture: no centroid dot, no connecting line -- just the
# polygon-inflated buffer contour + intersection region + dots colored
# by capture status.
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
# RENDER EACH PAIR
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

  # Recompute radius using our POLYGON-based find_capture_radius so
  # the plotted buffer follows the boundary of inner rather than a circle.
  ctarget <- find_capture_radius(inner, outer, res$dots,
                                 target_fraction = 0.50)

  plot_capture(
    inner, outer, res$dots,
    captured     = ctarget$captured,
    buffer       = ctarget$buffer,
    intersection = ctarget$intersection,
    title    = sprintf("Radius for 50%% capture"),
    subtitle = sprintf("r = %.2f km  |  achieved = %.3f",
                       ctarget$radius / 1000, ctarget$proportion)
  )
}

cat("\nRendering WWTP12 -> WWTP13 (close pair) ...\n")
p_left  <- render_pair_plot("WWTP12", "WWTP13")

cat("Rendering WWTP1 -> WWTP6 (moderate pair) ...\n")
p_right <- render_pair_plot("WWTP1", "WWTP6")

# ============================================================================
# COMBINE + SAVE
# ============================================================================
combined <- (p_left | p_right) +
  plot_annotation(
    caption = paste("Buffer expands from the boundary of the source polygon",
                    "(Omega_i) outward, following its irregular shape.")
  )

out_path <- file.path(HERE, "figure_2_3_combined_v2.pdf")
ggsave(out_path, combined, width = 12, height = 6)
cat("\n  -> ", out_path, "\n", sep = "")

cat("\nDone.\n")
