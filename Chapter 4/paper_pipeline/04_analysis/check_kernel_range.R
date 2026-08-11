# =============================================================================
# check_kernel_range.R
# -----------------------------------------------------------------------------
# Verifies the claim made in Section 2.3.2 of the manuscript that the Gaussian
# RBF kernel
#     K_ij = exp(-d_ij^2 / 2 sigma_K^2)
# applied to the PEHD matrix is bounded in the half-open interval (0, 1]:
#
#   - K_ii = 1 exactly (because d_ii = 0), so the upper bound 1 is achieved
#     on the diagonal => the closed bracket on the right is correct.
#   - All off-diagonal K_ij are strictly positive and strictly less than 1,
#     so the lower bound 0 is open and the upper bound is achieved only on
#     the diagonal.
#
# Run from the repo root:
#   Rscript "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/check_kernel_range.R"
# =============================================================================

PROJ_DIR  <- "/Users/josepalacio/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW"
PEHD_FILE <- file.path(PROJ_DIR, "Modelos AR con rho est. y rho=0",
                       "paper_pipeline", "results", "d_PEHD.rds")

d <- readRDS(PEHD_FILE)
P <- nrow(d)

# Paper's bandwidth formula: sigma_K^2 = 0.5 * median(d^2)
sigma_K_sq <- 0.5 * median(d[upper.tri(d)]^2)
sigma_K    <- sqrt(sigma_K_sq)
K          <- exp(-d^2 / (2 * sigma_K_sq))

diag_K  <- diag(K)
offdiag <- K[!diag(P)]

cat("============================================================\n")
cat("  Verifying K_ij range in (0, 1]\n")
cat("============================================================\n")
cat(sprintf("\nBandwidth (paper formula sigma_K^2 = 0.5 * median(d^2)):\n"))
cat(sprintf("  sigma_K^2 = %.4f km^2,   sigma_K = %.4f km\n",
            sigma_K_sq, sigma_K))

cat(sprintf("\nDIAGONAL (%d values):\n", P))
cat(sprintf("  K_ii min   = %.10f\n", min(diag_K)))
cat(sprintf("  K_ii max   = %.10f\n", max(diag_K)))
cat(sprintf("  all K_ii == 1 exactly?  %s\n", all(diag_K == 1)))

cat(sprintf("\nOFF-DIAGONAL (%d values):\n", length(offdiag)))
cat(sprintf("  K_ij min   = %.10f  (largest distance pair)\n", min(offdiag)))
cat(sprintf("  K_ij max   = %.10f  (smallest distance pair)\n", max(offdiag)))
cat(sprintf("  all K_ij > 0 strictly?  %s\n", all(offdiag > 0)))
cat(sprintf("  all K_ij < 1 strictly?  %s\n", all(offdiag < 1)))

cat(sprintf("\nGLOBAL:\n"))
cat(sprintf("  min(K) = %.10f   (> 0, so lower bound is OPEN)\n", min(K)))
cat(sprintf("  max(K) = %.10f   (= 1, only on diagonal => upper bound is CLOSED)\n",
            max(K)))
cat(sprintf("  any K_ij == 0?               %s\n", any(K == 0)))
cat(sprintf("  any off-diagonal K_ij == 1?  %s\n", any(offdiag == 1)))

cat("\n[Verdict]\n")
if (all(diag_K == 1) && all(offdiag > 0) && all(offdiag < 1)) {
  cat("  K_ij is bounded in (0, 1] : confirmed on the Houston PEHD data.\n")
  cat("  Lower bound OPEN  (no entry equals zero exactly).\n")
  cat("  Upper bound CLOSED (achieved on the diagonal: K_ii = 1).\n")
} else {
  cat("  Unexpected: claim NOT verified. Review the diagnostics above.\n")
}

cat("\n=== DONE ===\n")
