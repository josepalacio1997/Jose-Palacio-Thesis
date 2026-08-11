# =============================================================================
# check_pehd_negative_type.R
# -----------------------------------------------------------------------------
# Schoenberg negative-type test for the PEHD matrix.
#
# Theory:
#   The Gaussian RBF kernel K_ij = exp(-d_ij^2 / 2 sigma^2) is PSD for all
#   bandwidths sigma > 0 if and only if d is of negative type, equivalently
#   d^2 is conditionally negative definite (CND), equivalently the
#   double-centred Gram matrix
#       B = -1/2 * J * D2 * J,    J = I - (1/n) 11^T,    D2 = (d_ij^2)
#   is positive semi-definite. If B has any strictly negative eigenvalue,
#   d cannot be embedded isometrically into any L^2 space, and the
#   Schoenberg result does NOT guarantee K to be PSD --- not for any sigma.
#
#   The eigenvalues of B also give the classical MDS spectrum: the number
#   of negative eigenvalues is the embedding dimension shortfall, and
#   their magnitude quantifies how far d is from being Euclidean.
#
# Run from the repo root:
#   Rscript "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/check_pehd_negative_type.R"
# =============================================================================

PROJ_DIR  <- "/Users/josepalacio/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW"
PEHD_FILE <- file.path(PROJ_DIR, "Modelos AR con rho est. y rho=0",
                       "paper_pipeline", "results", "d_PEHD.rds")
OUT_CSV   <- file.path(PROJ_DIR, "Modelos AR con rho est. y rho=0",
                       "paper_pipeline", "results",
                       "pehd_negative_type_B_eigenvalues.csv")
WRITE_CSV <- TRUE

d <- readRDS(PEHD_FILE)
P <- nrow(d)
D2 <- d^2

cat("============================================================\n")
cat("  Schoenberg negative-type test on PEHD\n")
cat("============================================================\n\n")
cat(sprintf("PEHD matrix: %d x %d\n", P, P))
cat(sprintf("  d range:    [%.3f, %.3f] km\n",
            min(d[upper.tri(d)]), max(d[upper.tri(d)])))
cat(sprintf("  d^2 range:  [%.3f, %.3f] km^2\n\n",
            min(D2[upper.tri(D2)]), max(D2[upper.tri(D2)])))

# ============================================================================
# 1. Double-centred Gram matrix B = -1/2 * J * D2 * J
# ============================================================================
J <- diag(P) - matrix(1 / P, P, P)
B <- -0.5 * J %*% D2 %*% J
B <- 0.5 * (B + t(B))            # enforce numerical symmetry

eig_B <- eigen(B, symmetric = TRUE, only.values = TRUE)$values

cat("[Double-centred Gram matrix B = -1/2 J D^(2) J]\n")
cat(sprintf("  B is %d x %d, trace(B) = %.4f\n", P, P, sum(diag(B))))
cat(sprintf("  eigenvalue range: [%.4f, %.4f]\n", min(eig_B), max(eig_B)))
cat(sprintf("  negative eigenvalues: %d / %d\n",
            sum(eig_B < -1e-10), length(eig_B)))
cat(sprintf("  sum of negatives: %+.4f   sum of positives: %+.4f\n",
            sum(eig_B[eig_B < 0]), sum(eig_B[eig_B > 0])))
cat("  full spectrum (sorted ascending):\n")
print(round(sort(eig_B), 6))
cat("\n")

# ============================================================================
# 2. Verdict
# ============================================================================
is_psd <- all(eig_B > -1e-10)
cat("[Schoenberg verdict]\n")
if (is_psd) {
  cat("  B is PSD  =>  d is of NEGATIVE TYPE  =>  K is PSD for all sigma > 0.\n")
  cat("  d embeds isometrically into R^", sum(eig_B > 1e-10), "\n", sep = "")
} else {
  shortfall <- sum(eig_B < -1e-10)
  cat(sprintf("  B has %d strictly negative eigenvalues.\n", shortfall))
  cat("  => d is NOT of negative type.\n")
  cat("  => d does NOT embed isometrically into any Hilbert space.\n")
  cat("  => Schoenberg's theorem does NOT guarantee K to be PSD for any sigma.\n")
  cat(sprintf("  Magnitude of non-Euclideanness: |sum of B's negative eigs| = %.4f km^2,\n",
              abs(sum(eig_B[eig_B < 0]))))
  cat(sprintf("    which is %.2f%% of |sum of B's positive eigs|.\n",
              100 * abs(sum(eig_B[eig_B < 0])) / sum(eig_B[eig_B > 0])))
}
cat("\n")

# ============================================================================
# 3. (Optional) Save eigenvalues of B for the appendix
# ============================================================================
if (WRITE_CSV) {
  out <- data.frame(rank = seq_along(eig_B), eigenvalue = sort(eig_B))
  write.csv(out, OUT_CSV, row.names = FALSE)
  cat(sprintf("[CSV] wrote: %s\n", OUT_CSV))
}

cat("\n=== DONE ===\n")
