# ============================================================================
# compute_W_matrices.R
# -----------------------------------------------------------------------------
# Computes the THREE spatial weight matrices used in the paper:
#   - W_PEHD  : Population-Extended Hausdorff Distance (paper's primary)
#   - W_eucl  : Euclidean centroid-to-centroid distance (ablation #1)
#   - W_HD    : Classical Hausdorff distance (ablation #2)
#
# All three go through the SAME downstream pipeline:
#   distance d_ij -> RBF kernel K_ij -> zero diagonal -> row-normalize -> W_ij
#
# Median heuristic bandwidth is computed independently for each distance
# matrix, so the comparison is fair: each kernel uses its own natural scale.
#
# Outputs (saved to ../results/):
#   W_PEHD.rds, d_PEHD.rds   : 15 x 15 matrices (rows of W sum to 1)
#   W_eucl.rds, d_eucl.rds   : ditto
#   W_HD.rds,   d_HD.rds     : ditto
# ============================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(readr); library(tibble)
})

# Project paths (this script lives in paper_pipeline/01_W_matrix/)
# Run from this directory: Rscript compute_W_matrices.R
PIPE_ROOT   <- "/Users/josepalacio/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW"
RESULTS_DIR <- file.path(PIPE_ROOT, "Modelos AR con rho est. y rho=0",
                         "paper_pipeline", "results")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# Source the distance machinery
source(file.path(PIPE_ROOT, "R", "core.R"))

# ---- Inputs -----------------------------------------------------------------
shp_batch       <- file.path(PIPE_ROOT, "data", "WWTP_Batch",
                             "WWTP_061621_mergedWWTPs to Batch_updatedon07192024.shp")
tracts_shp_path <- file.path(PIPE_ROOT, "data", "Census_Tracts_Population",
                             "WWTP Batch Service Areas_Census Tracts_Population.shp")

Batch_focus <- paste0("WWTP", 1:15)

# ---- Load and prepare polygons ----------------------------------------------
wwtp_shp_batch <- sf::read_sf(shp_batch) |>
  sf::st_transform(32615) |>
  mutate(BatchName = toupper(trimws(BatchName)))

batch_key_shp <- wwtp_shp_batch |>
  sf::st_drop_geometry() |> distinct(BatchName) |> arrange(BatchName) |>
  mutate(Batch = paste0("WWTP", row_number()))

wwtp_shp_batch <- wwtp_shp_batch |> left_join(batch_key_shp, by = "BatchName")
wwtp_shp_focus <- wwtp_shp_batch |>
  filter(Batch %in% Batch_focus) |>
  arrange(match(Batch, Batch_focus))

tracts_raw <- sf::st_read(tracts_shp_path, quiet = TRUE) |>
  sf::st_make_valid() |>
  dplyr::filter(sf::st_is_valid(geometry), !sf::st_is_empty(geometry))
tracts_valid <- tracts_raw |>
  dplyr::distinct(GEOID, .keep_all = TRUE) |>
  sf::st_transform(32615)

# ---- Helper: RBF + row-normalization ----------------------------------------
make_W_from_distance <- function(d, sigma = NULL) {
  d_dim   <- dim(d); d_names <- dimnames(d)
  d_num   <- structure(as.numeric(unclass(d)), dim = d_dim, dimnames = d_names)
  off     <- d_num; diag(off) <- NA
  if (is.null(sigma)) sigma <- sqrt(median(off^2, na.rm = TRUE) / 2)
  K <- exp(-d_num^2 / (2 * sigma^2)); diag(K) <- 0
  W <- K / rowSums(K)
  attr(W, "sigma_K") <- sigma
  W
}

bn_to_wwtp <- setNames(batch_key_shp$Batch, batch_key_shp$BatchName)

# ============================================================================
# 1. PEHD matrix
# ============================================================================
cat("\n=== PEHD (Population-Extended Hausdorff Distance) ===\n")
d_PEHD <- region_distance_matrix(
  regions         = wwtp_shp_focus,
  label_col       = "BatchName",
  method          = "population",
  tracts          = tracts_valid,
  pop_col         = "POP24_5YR",
  target_fraction = 0.5,
  people_per_dot  = 200,
  units           = "km"
)
W_PEHD <- make_W_from_distance(d_PEHD)
rownames(W_PEHD) <- bn_to_wwtp[rownames(W_PEHD)]
colnames(W_PEHD) <- bn_to_wwtp[colnames(W_PEHD)]
W_PEHD <- W_PEHD[Batch_focus, Batch_focus, drop = FALSE]

# Also re-name d_PEHD to WWTP labels for consistency
rownames(d_PEHD) <- bn_to_wwtp[rownames(d_PEHD)]
colnames(d_PEHD) <- bn_to_wwtp[colnames(d_PEHD)]
d_PEHD <- d_PEHD[Batch_focus, Batch_focus, drop = FALSE]

cat(sprintf("  sigma_K (PEHD)  : %.3f km\n", attr(W_PEHD, "sigma_K")))
cat(sprintf("  off-diag d range: [%.2f, %.2f] km\n",
            min(d_PEHD[lower.tri(d_PEHD)]), max(d_PEHD[lower.tri(d_PEHD)])))

# ============================================================================
# 2. Euclidean matrix (ablation)
# ============================================================================
cat("\n=== Euclidean (centroid-to-centroid) ===\n")
d_eucl <- region_distance_matrix(
  regions   = wwtp_shp_focus,
  label_col = "BatchName",
  method    = "centroid",
  units     = "km"
)
W_eucl <- make_W_from_distance(d_eucl)
rownames(W_eucl) <- bn_to_wwtp[rownames(W_eucl)]
colnames(W_eucl) <- bn_to_wwtp[colnames(W_eucl)]
W_eucl <- W_eucl[Batch_focus, Batch_focus, drop = FALSE]

rownames(d_eucl) <- bn_to_wwtp[rownames(d_eucl)]
colnames(d_eucl) <- bn_to_wwtp[colnames(d_eucl)]
d_eucl <- d_eucl[Batch_focus, Batch_focus, drop = FALSE]

cat(sprintf("  sigma_K (eucl)  : %.3f km\n", attr(W_eucl, "sigma_K")))
cat(sprintf("  off-diag d range: [%.2f, %.2f] km\n",
            min(d_eucl[lower.tri(d_eucl)]), max(d_eucl[lower.tri(d_eucl)])))

# ============================================================================
# 3. Hausdorff matrix (classical Hausdorff distance, ablation #2)
# ============================================================================
cat("\n=== HD (Classical Hausdorff Distance) ===\n")
d_HD <- region_distance_matrix(
  regions   = wwtp_shp_focus,
  label_col = "BatchName",
  method    = "hausdorff",
  units     = "km"
)
W_HD <- make_W_from_distance(d_HD)
rownames(W_HD) <- bn_to_wwtp[rownames(W_HD)]
colnames(W_HD) <- bn_to_wwtp[colnames(W_HD)]
W_HD <- W_HD[Batch_focus, Batch_focus, drop = FALSE]

rownames(d_HD) <- bn_to_wwtp[rownames(d_HD)]
colnames(d_HD) <- bn_to_wwtp[colnames(d_HD)]
d_HD <- d_HD[Batch_focus, Batch_focus, drop = FALSE]

cat(sprintf("  sigma_K (HD)    : %.3f km\n", attr(W_HD, "sigma_K")))
cat(sprintf("  off-diag d range: [%.2f, %.2f] km\n",
            min(d_HD[lower.tri(d_HD)]), max(d_HD[lower.tri(d_HD)])))

# ============================================================================
# 4. Save
# ============================================================================
saveRDS(W_PEHD, file.path(RESULTS_DIR, "W_PEHD.rds"))
saveRDS(d_PEHD, file.path(RESULTS_DIR, "d_PEHD.rds"))
saveRDS(W_eucl, file.path(RESULTS_DIR, "W_eucl.rds"))
saveRDS(d_eucl, file.path(RESULTS_DIR, "d_eucl.rds"))
saveRDS(W_HD,   file.path(RESULTS_DIR, "W_HD.rds"))
saveRDS(d_HD,   file.path(RESULTS_DIR, "d_HD.rds"))

cat("\n=== DONE ===\n")
cat(sprintf("Wrote: %s/{W_PEHD, d_PEHD, W_eucl, d_eucl, W_HD, d_HD}.rds\n",
            RESULTS_DIR))

# ============================================================================
# 5. Quick numerical comparison
# ============================================================================
cat("\n=== Pairwise comparison of the three W matrices ===\n")
cat(sprintf("cor(vec(W_PEHD), vec(W_eucl))  : %.4f\n",
            cor(as.numeric(W_PEHD), as.numeric(W_eucl))))
cat(sprintf("cor(vec(W_PEHD), vec(W_HD))    : %.4f\n",
            cor(as.numeric(W_PEHD), as.numeric(W_HD))))
cat(sprintf("cor(vec(W_eucl), vec(W_HD))    : %.4f\n",
            cor(as.numeric(W_eucl), as.numeric(W_HD))))
cat(sprintf("max |W_PEHD - W_eucl|          : %.4f\n",
            max(abs(W_PEHD - W_eucl))))
cat(sprintf("max |W_PEHD - W_HD|            : %.4f\n",
            max(abs(W_PEHD - W_HD))))
cat(sprintf("max |W_eucl - W_HD|            : %.4f\n",
            max(abs(W_eucl - W_HD))))
