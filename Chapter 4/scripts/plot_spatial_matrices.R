# plot_spatial_matrices.R --------------------------------------------------
# Reproduces Figures 2 and 3 of Palacio et al. (2026) using the output
# of build_spatial_matrices().
#
# Outputs (in Modelos AR con rho est. y rho=0/):
#   figure_2_pehd_example.pdf        — asymmetric PEHD illustration
#                                       (WWTP12 -> WWTP13, WWTP1 -> WWTP3)
#   figure_3_D_heatmap.pdf           — left panel of paper Figure 3
#   figure_3_W_heatmap.pdf           — right panel of paper Figure 3
#   figure_3_combined.pdf            — both heatmaps side by side
#
# Run:  Rscript scripts/plot_spatial_matrices.R
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(patchwork)
})

# -------------------------------------------------------------------- Paths
PROJ_ROOT <- if (Sys.info()["sysname"] == "Darwin") {
  "~/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW"
} else {
  normalizePath(".", mustWork = FALSE)
}
PROJ_ROOT <- normalizePath(PROJ_ROOT, mustWork = TRUE)

OUT_DIR <- file.path(PROJ_ROOT, "Modelos AR con rho est. y rho=0",
                     "figures_from_build_spatial")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

shp_batch       <- file.path(PROJ_ROOT, "data", "WWTP_Batch",
  "WWTP_061621_mergedWWTPs to Batch_updatedon07192024.shp")
tracts_shp_path <- file.path(PROJ_ROOT, "data", "Census_Tracts_Population",
  "WWTP Batch Service Areas_Census Tracts_Population.shp")

stopifnot(file.exists(shp_batch), file.exists(tracts_shp_path))

# ---------------------------------------------------------------- Load data
cat("Loading data...\n")
Batch_focus <- paste0("WWTP", 1:15)

wwtp_shp <- sf::read_sf(shp_batch) |>
  sf::st_transform(32615) |>
  mutate(BatchName = toupper(trimws(BatchName)))

batch_key <- wwtp_shp |>
  sf::st_drop_geometry() |>
  distinct(BatchName) |>
  arrange(BatchName) |>
  mutate(Batch = paste0("WWTP", row_number()))

wwtp_shp_focus <- wwtp_shp |>
  left_join(batch_key, by = "BatchName") |>
  filter(Batch %in% Batch_focus) |>
  arrange(match(Batch, Batch_focus))

tracts_raw <- sf::st_read(tracts_shp_path, quiet = TRUE) |>
  sf::st_make_valid() |>
  dplyr::filter(sf::st_is_valid(geometry), !sf::st_is_empty(geometry))

tracts_valid <- tracts_raw |>
  dplyr::distinct(GEOID, .keep_all = TRUE) |>
  sf::st_transform(32615)

# --------------------------------------- Run build_spatial_matrices()
cat("Running build_spatial_matrices()...\n")
source(file.path(PROJ_ROOT, "R", "build_spatial_matrices.R"))
spatial <- build_spatial_matrices(
  regions         = wwtp_shp_focus,
  tracts          = tracts_valid,
  pop_col         = "POP24_5YR",
  people_per_dot  = 200,
  target_fraction = 0.5,
  labels          = Batch_focus
)

# ============================================================================
# FIGURE 3 — Heatmaps of D and W
# ============================================================================

cat("\nGenerating Figure 3 (D and W heatmaps)...\n")

# Convert matrices to long-form data frames for ggplot
mat_to_long <- function(M, value_name = "value") {
  df <- as.data.frame(as.table(M))
  names(df) <- c("source", "target", value_name)
  df$source <- factor(df$source, levels = Batch_focus)
  df$target <- factor(df$target, levels = Batch_focus)
  df
}

df_D <- mat_to_long(spatial$D, "km")
df_W <- mat_to_long(spatial$W, "Wij")

# Common plotting theme (matches paper style)
heatmap_theme <- theme_minimal(base_size = 10) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y        = element_text(size = 8),
    axis.title         = element_blank(),
    panel.grid         = element_blank(),
    legend.position    = "bottom",
    legend.key.height  = unit(0.3, "cm"),
    legend.key.width   = unit(0.9, "cm"),
    plot.title         = element_text(hjust = 0.5, face = "bold")
  )

# Left panel: D (PEHD in km) — WWTP1 at top-left
p_D <- ggplot(df_D, aes(x = target, y = source, fill = km)) +
  geom_tile(color = "white", size = 0.15) +
  scale_fill_viridis_c(name = "km", option = "magma", direction = -1,
                       na.value = "grey90") +
  scale_y_discrete(limits = rev(Batch_focus)) +
  coord_equal() +
  labs(title = "D  (PEHD, km)") +
  heatmap_theme

# Right panel: W (row-normalized weights) — WWTP1 at top-left
p_W <- ggplot(df_W, aes(x = target, y = source, fill = Wij)) +
  geom_tile(color = "white", size = 0.15) +
  scale_fill_viridis_c(name = expression(W[ij]), option = "magma",
                       direction = -1, na.value = "grey90",
                       limits = c(0, max(df_W$Wij))) +
  scale_y_discrete(limits = rev(Batch_focus)) +
  coord_equal() +
  labs(title = "W  (spatial weights)") +
  heatmap_theme

# Save individual + combined
ggsave(file.path(OUT_DIR, "figure_3_D_heatmap.pdf"),
       p_D, width = 6, height = 6)
ggsave(file.path(OUT_DIR, "figure_3_W_heatmap.pdf"),
       p_W, width = 6, height = 6)
ggsave(file.path(OUT_DIR, "figure_3_combined.pdf"),
       p_D + p_W, width = 12, height = 6.5)
cat("  → figure_3_D_heatmap.pdf\n")
cat("  → figure_3_W_heatmap.pdf\n")
cat("  → figure_3_combined.pdf\n")

# ============================================================================
# FIGURE 2 — Asymmetric PEHD illustration
# ============================================================================

cat("\nGenerating Figure 2 (asymmetric PEHD examples)...\n")

# Helper: compute the achieved fraction at a given radius r (source Omega_i
# is regions[i], target dots are dots_per_region[[j]]).
achieved_fraction <- function(i, j, r_km) {
  dots  <- spatial$dots_per_region[[j]]
  dists <- as.numeric(sf::st_distance(wwtp_shp_focus[i, ], dots)) / 1000
  mean(dists <= r_km)
}

# Build one panel for a source-target pair.
make_pehd_panel <- function(i, j, title_prefix = "") {
  src   <- wwtp_shp_focus[i, ]
  tgt   <- wwtp_shp_focus[j, ]
  dots  <- spatial$dots_per_region[[j]]

  r_km  <- spatial$D_tilde[i, j]         # asymmetric distance i -> j
  frac  <- achieved_fraction(i, j, r_km)

  # Expand a buffer of r_km around the source polygon (in meters)
  buf <- sf::st_buffer(src, dist = r_km * 1000)

  # Classify dots as inside/outside the buffer
  dots_sf <- sf::st_as_sf(dots)
  inside  <- lengths(sf::st_intersects(dots_sf, buf)) > 0
  dots_sf$status <- factor(ifelse(inside, "inside", "outside"),
                           levels = c("inside", "outside"))

  # Transform back to lat/lon for a friendlier map view
  src4326 <- sf::st_transform(src,     4326)
  tgt4326 <- sf::st_transform(tgt,     4326)
  buf4326 <- sf::st_transform(buf,     4326)
  dts4326 <- sf::st_transform(dots_sf, 4326)

  title <- sprintf("%s%s → %s\n r = %.2f km | achieved = %.3f",
                   title_prefix,
                   Batch_focus[i], Batch_focus[j], r_km, frac)

  ggplot() +
    geom_sf(data = tgt4326,             fill = "lightblue", color = "steelblue",
            size = 0.4, alpha = 0.5) +
    geom_sf(data = src4326,             fill = "salmon",    color = "firebrick",
            size = 0.4, alpha = 0.5) +
    geom_sf(data = buf4326, fill = NA, color = "black", linetype = "dashed",
            linewidth = 0.5) +
    geom_sf(data = dts4326, aes(color = status), size = 0.25) +
    scale_color_manual(values = c(inside = "gold", outside = "steelblue4"),
                       name = "Dot status") +
    labs(title = title) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 9),
          legend.position = "bottom",
          axis.text = element_text(size = 7))
}

# Paper Figure 2 uses (WWTP12 → WWTP13) [close pair] and
# (WWTP1 → WWTP3) [distant pair].
p_close  <- make_pehd_panel(i = 12, j = 13)
p_distant <- make_pehd_panel(i = 1,  j = 3)

fig2 <- p_close + p_distant + plot_layout(ncol = 2)
ggsave(file.path(OUT_DIR, "figure_2_pehd_example.pdf"),
       fig2, width = 12, height = 6.5)
cat("  → figure_2_pehd_example.pdf\n")

cat("\nAll figures written to:\n  ", OUT_DIR, "\n", sep = "")
cat("Done.\n")
