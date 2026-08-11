# =============================================================================
# check_pehd_diagnostics.R
# -----------------------------------------------------------------------------
# Standalone PEHD diagnostic script. Reads the PEHD distance matrix produced by
# 01_W_matrix/compute_W_matrices.R and verifies, on the actual Houston data:
#
#   1. Basic shape: symmetry, zero diagonal, off-diagonal range.
#   2. Triangle inequality: enumerate all ordered triples (i, j, k) and
#      count how many violate d_{ik} <= d_{ij} + d_{jk}. Reports the worst
#      offenders.
#   3. Kernel spectrum: build K_{ij} = exp(-d_{ij}^2 / 2 sigma_K^2) using the
#      paper's bandwidth (sigma_K^2 = 0.5 * median{d_{ij}^2 : i != j}) and
#      report the full eigenvalue spectrum, range, and count of negatives.
#
# Optionally writes:
#   paper_pipeline/results/pehd_triangle_violations.csv  (all violating triples)
#   paper_pipeline/results/kernel_K_eigenvalues.csv      (15 sorted eigenvalues)
#
# Run from the repo root:
#   Rscript "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/check_pehd_diagnostics.R"
# =============================================================================

PROJ_DIR <- "/Users/josepalacio/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW"
PEHD_FILE <- file.path(PROJ_DIR, "Modelos AR con rho est. y rho=0",
                       "paper_pipeline", "results", "d_PEHD.rds")
RESULTS_DIR <- file.path(PROJ_DIR, "Modelos AR con rho est. y rho=0",
                         "paper_pipeline", "results")
WRITE_CSV <- TRUE   # set FALSE to skip CSV output

d <- readRDS(PEHD_FILE)
P <- nrow(d)

cat("============================================================\n")
cat("  PEHD diagnostic: triangle inequality + kernel spectrum\n")
cat("============================================================\n\n")

# ============================================================================
# 1. Basic shape
# ============================================================================
cat(sprintf("d matrix: %d x %d\n", P, P))
cat(sprintf("  symmetry max |d - t(d)|: %.6g\n", max(abs(d - t(d)))))
cat(sprintf("  diagonal max |d_ii|:     %.6g\n", max(abs(diag(d)))))
cat(sprintf("  off-diagonal range:      [%.3f, %.3f] km\n\n",
            min(d[upper.tri(d)]), max(d[upper.tri(d)])))

# ============================================================================
# 2. Triangle inequality scan
# ============================================================================
viol <- list()
for (i in 1:P) for (j in 1:P) for (k in 1:P) {
  if (i == j || j == k || i == k) next
  excess <- d[i, k] - (d[i, j] + d[j, k])
  if (excess > 1e-9) {
    viol[[length(viol) + 1]] <- data.frame(
      i = i, j = j, k = k,
      d_ik         = d[i, k],
      d_ij_plus_jk = d[i, j] + d[j, k],
      excess       = excess
    )
  }
}
viol_df <- do.call(rbind, viol)
total   <- P * (P - 1) * (P - 2)

cat("[Triangle inequality]\n")
cat(sprintf("  violating triples: %d / %d  (%.1f%%)\n",
            nrow(viol_df), total, 100 * nrow(viol_df) / total))
cat(sprintf("  worst excess:      %.3f km\n", max(viol_df$excess)))
cat("  top 5 worst triples (i, j, k):\n")
print(head(viol_df[order(-viol_df$excess), ], 5), row.names = FALSE)
cat("\n")

# ============================================================================
# 3. Kernel spectrum (using the paper's exact bandwidth formula)
# ============================================================================
sigma_K_sq <- 0.5 * median(d[upper.tri(d)]^2)
sigma_K    <- sqrt(sigma_K_sq)
K   <- exp(-d^2 / (2 * sigma_K_sq))
eig <- eigen(K, symmetric = TRUE, only.values = TRUE)$values

cat("[Gaussian RBF kernel, paper formula sigma_K^2 = 0.5 * median(d^2)]\n")
cat(sprintf("  sigma_K^2 = %.4f km^2,   sigma_K = %.4f km\n",
            sigma_K_sq, sigma_K))
cat(sprintf("  K is %d x %d, trace(K) = %.4f\n", P, P, sum(diag(K))))
cat(sprintf("  eigenvalue range: [%.4f, %.4f]\n", min(eig), max(eig)))
cat(sprintf("  negative eigenvalues: %d / %d\n",
            sum(eig < 0), length(eig)))
cat(sprintf("  sum of negatives: %+.4f   sum of positives: %+.4f\n",
            sum(eig[eig < 0]), sum(eig[eig > 0])))
cat("  full spectrum (sorted):\n")
print(round(sort(eig), 6))
cat("\n")

# ============================================================================
# 4. Optional: write CSVs for paper / appendix
# ============================================================================
if (WRITE_CSV) {
  out_viol <- file.path(RESULTS_DIR, "pehd_triangle_violations.csv")
  out_eig  <- file.path(RESULTS_DIR, "kernel_K_eigenvalues.csv")
  write.csv(viol_df[order(-viol_df$excess), ], out_viol, row.names = FALSE)
  write.csv(
    data.frame(rank = seq_along(eig), eigenvalue = sort(eig)),
    out_eig, row.names = FALSE
  )
  cat(sprintf("[CSV] wrote: %s\n", out_viol))
  cat(sprintf("[CSV] wrote: %s\n", out_eig))
}

cat("\n=== DONE ===\n")
