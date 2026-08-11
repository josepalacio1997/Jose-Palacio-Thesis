# test_build_spatial_matrices.R --------------------------------------------
# Quick smoke-test of build_spatial_matrices() on the real Houston data.
# Verifies:
#   1. Function runs end-to-end without errors.
#   2. Output matrices have the expected shape (15 x 15).
#   3. Sanity properties hold: rowSums(W) == 1, diag(W) == 0, D symmetric.
#   4. Numbers match those in the paper (Figure 3 caption ranges).
#
# Run:  Rscript scripts/test_build_spatial_matrices.R
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
})

# -------------------------------------------------------------------- Paths
PROJ_ROOT <- if (Sys.info()["sysname"] == "Darwin") {
  "~/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW"
} else {
  normalizePath(".", mustWork = FALSE)
}
PROJ_ROOT <- normalizePath(PROJ_ROOT, mustWork = TRUE)

shp_batch       <- file.path(PROJ_ROOT, "data", "WWTP_Batch",
  "WWTP_061621_mergedWWTPs to Batch_updatedon07192024.shp")
tracts_shp_path <- file.path(PROJ_ROOT, "data", "Census_Tracts_Population",
  "WWTP Batch Service Areas_Census Tracts_Population.shp")

stopifnot(file.exists(shp_batch), file.exists(tracts_shp_path))

# ---------------------------------------------------------------- Load WWTPs
cat("Loading WWTP polygons...\n")
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

stopifnot(nrow(wwtp_shp_focus) == 15)
cat(sprintf("  %d WWTP polygons loaded.\n", nrow(wwtp_shp_focus)))

# -------------------------------------------------------- Load census tracts
cat("Loading census tracts...\n")
tracts_raw <- sf::st_read(tracts_shp_path, quiet = TRUE) |>
  sf::st_make_valid() |>
  dplyr::filter(sf::st_is_valid(geometry), !sf::st_is_empty(geometry))

tracts_valid <- tracts_raw |>
  dplyr::distinct(GEOID, .keep_all = TRUE) |>
  sf::st_transform(32615)

cat(sprintf("  %d census tracts (deduped, valid geometry).\n",
            nrow(tracts_valid)))

# -------------------------------------- Source the function under test
cat("\nSourcing build_spatial_matrices()...\n")
source(file.path(PROJ_ROOT, "R", "build_spatial_matrices.R"))

# ----------------------------------------------------------- Run the function
cat("Running build_spatial_matrices() ...\n")
t0 <- Sys.time()

spatial <- build_spatial_matrices(
  regions         = wwtp_shp_focus,
  tracts          = tracts_valid,
  pop_col         = "POP24_5YR",
  people_per_dot  = 200,
  target_fraction = 0.5,
  labels          = Batch_focus
)

elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  build_spatial_matrices() completed in %.1f seconds.\n",
            elapsed))

# ---------------------------------------------------------- Basic shape check
cat("\n=== Shape check ===\n")
for (nm in c("D_tilde", "D", "K", "W")) {
  M <- spatial[[nm]]
  cat(sprintf("  %-8s: %d x %d, dimnames(1) = %s\n",
              nm, nrow(M), ncol(M), rownames(M)[1]))
  stopifnot(dim(M) == c(15, 15))
}

# ----------------------------------------------------- Sanity property checks
cat("\n=== Property checks ===\n")

cat(sprintf("  D symmetric:                        %s\n",
            isSymmetric(unname(spatial$D))))
cat(sprintf("  diag(D) all zero:                   %s\n",
            all(diag(spatial$D) == 0)))
cat(sprintf("  D >= D_tilde and >= t(D_tilde):     %s\n",
            all(spatial$D >= spatial$D_tilde) &&
            all(spatial$D >= t(spatial$D_tilde))))
cat(sprintf("  K symmetric (before zero diag):     %s (informational)\n",
            "K has zero diag now"))
cat(sprintf("  diag(K) all zero:                   %s\n",
            all(diag(spatial$K) == 0)))
cat(sprintf("  rowSums(W) == 1 (within 1e-10):     %s\n",
            all(abs(rowSums(spatial$W) - 1) < 1e-10)))
cat(sprintf("  diag(W) all zero:                   %s\n",
            all(diag(spatial$W) == 0)))

# -------------------------------------------- Numerical values (paper Fig. 3)
cat("\n=== Numerical values (compare to paper Figure 3) ===\n")
off_D  <- spatial$D[row(spatial$D) != col(spatial$D)]
cat(sprintf("  D off-diagonal range: [%.2f, %.2f] km\n",
            min(off_D), max(off_D)))
cat(sprintf("     Paper Figure 3 says: 'approximately 2 km to ~49 km'\n"))

cat(sprintf("  sigma_K (median heuristic): %.2f km\n", spatial$sigma_K))
off_W <- spatial$W[row(spatial$W) != col(spatial$W)]
cat(sprintf("  W off-diagonal range: [%.4f, %.4f]\n",
            min(off_W), max(off_W)))
cat(sprintf("     Paper Figure 3 legend: 'W_ij in [0.0, 0.4]'\n"))

# ------------------------------------------------------- Dot counts per plant
cat("\n=== Dot counts per WWTP (M_j) ===\n")
dot_counts <- vapply(spatial$dots_per_region, length, 0L)
names(dot_counts) <- Batch_focus
print(dot_counts)
cat(sprintf("  Total dots: %d (people_per_dot = %d)\n",
            sum(dot_counts), spatial$people_per_dot))

# --------------------------------------------------- Compare with legacy code
cat("\n=== Sanity vs legacy `R/core.R` (optional) ===\n")
if (file.exists(file.path(PROJ_ROOT, "R", "core.R"))) {
  source(file.path(PROJ_ROOT, "R", "core.R"))
  cat("  Running legacy region_distance_matrix() for cross-check...\n")
  d_legacy <- region_distance_matrix(
    regions = wwtp_shp_focus, label_col = "BatchName",
    method = "population", tracts = tracts_valid,
    pop_col = "POP24_5YR", target_fraction = 0.5,
    people_per_dot = 200, units = "km"
  )
  # legacy has BatchName rownames; align by reordering
  bn_to_wwtp <- setNames(batch_key$Batch, batch_key$BatchName)
  rownames(d_legacy) <- bn_to_wwtp[rownames(d_legacy)]
  colnames(d_legacy) <- bn_to_wwtp[colnames(d_legacy)]
  d_legacy <- d_legacy[Batch_focus, Batch_focus]

  d_new    <- spatial$D
  max_diff <- max(abs(d_new - d_legacy))
  rel_diff <- max_diff / max(abs(d_legacy))
  cat(sprintf("  Max absolute diff  new D vs legacy: %.4e km\n", max_diff))
  cat(sprintf("  Max relative diff  new D vs legacy: %.4e\n",   rel_diff))
  if (rel_diff < 0.05) {
    cat("  New function agrees with legacy within 5% (sampling stochastic).\n")
  } else {
    cat("  WARNING: relative diff exceeds 5%. Investigate.\n")
  }
} else {
  cat("  R/core.R not found; skipping legacy cross-check.\n")
}

cat("\n=== Done. ===\n")
