# =============================================================================
# check_sigma_obs_partition.R
# -----------------------------------------------------------------------------
# Empirical check of the identifiability claim in the Discussion: does the
# Direct variant attribute LESS sigma_obs than the SSM variant, pushing more
# variance into the spatial coupling parameter rho?
#
# To compare apples-to-apples we need BOTH:
#   - SSM:    sigma_filt(p, t) from the MARSS Stage-1 smoother (data passed
#             to Stan, not saved in the fit). We rebuild it from y_batch
#             using the EXACT spec from run_A_marss_est.R.
#   - Direct: sigma_obs_pl[p] posterior median (parameter estimated in Stan;
#             extracted from the saved draws).
#
# Output: side-by-side per-plant table + aggregate ratio + CSV artifact.
#
# Run from the repo root:
#   Rscript "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/check_sigma_obs_partition.R"
#
# Runtime: ~30-60 seconds (15 MARSS fits + draws subset).
# =============================================================================

suppressPackageStartupMessages({
  library(MARSS); library(posterior); library(dplyr); library(tidyr)
})

PROJ_DIR <- "/Users/josepalacio/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW"
FITS_DIR <- file.path(PROJ_DIR, "Modelos AR con rho est. y rho=0", "fits")

P <- 15
plant_names <- paste0("WWTP", 1:P)
ln10 <- log(10)

# ============================================================================
# 1 . Load y_raw matrix from any fit (P x T, natural scale = B gc/day).
# ============================================================================
direct <- readRDS(file.path(FITS_DIR, "A_direct_est.rds"))
y_PT   <- direct$y_raw
if (is.null(y_PT)) stop("Could not find y_raw inside A_direct_est.rds.")
T_focus <- ncol(y_PT)
cat(sprintf("y_PT shape: %d plants x %d weeks\n", P, T_focus))

# Convert to log10 in T x P orientation (matching runner's Y_log10_TP).
Y_log10_TP <- log10(t(y_PT))
Y_log10_TP[!is.finite(Y_log10_TP)] <- NA

# ============================================================================
# 2 . MARSS Stage 1 (same smoother spec as run_A_marss_est.R: 2-state with
#     B = [[2,-1],[1,0]], R diagonal, Z = (1,0)).
# ============================================================================
ssm_spec_one_plant <- list(
  B  = matrix(c(2, -1, 1, 0), 2, 2, byrow = TRUE),
  U  = "zero", Q = "diagonal and equal", V0 = "identity", x0 = "equal",
  R  = "diagonal and equal", A = "zero", Z = matrix(c(1, 0), 1, 2)
)

sigma_filt <- matrix(NA_real_, P, T_focus,
                     dimnames = list(plant_names, NULL))

cat("Running MARSS Stage 1 (15 plant-level smoothers)...\n")
for (p in 1:P) {
  yp_log10 <- Y_log10_TP[, p]
  if (sum(!is.na(yp_log10)) < 5) {
    sigma_filt[p, ] <- 1.0
    next
  }
  fit_p <- tryCatch(
    MARSS::MARSS(matrix(yp_log10, nrow = 1),
                 model   = ssm_spec_one_plant,
                 method  = "BFGS", silent = TRUE,
                 control = list(maxit = 200)),
    error = function(e) NULL
  )
  if (is.null(fit_p)) {
    cat(sprintf("  WWTP%d: MARSS failed; sigma_filt set to 1.0\n", p))
    sigma_filt[p, ] <- 1.0
    next
  }
  kfas_p <- MARSS::MARSSkfas(fit_p)
  Vtt_p  <- kfas_p$Vtt                                    # 2 x 2 x T
  Mp <- coef(fit_p, type = "matrix")
  Zp <- as.matrix(Mp$Z)                                    # 1 x 2

  V_obs_t <- vapply(seq_len(T_focus),
                    function(tt) as.numeric(Zp %*% Vtt_p[, , tt] %*% t(Zp)),
                    numeric(1))
  V_obs_t  <- pmax(V_obs_t, 0)
  var_ln_t <- (ln10^2) * V_obs_t
  sigma_filt[p, ] <- sqrt(var_ln_t)
  cat(sprintf("  WWTP%-3d ok. sigma_filt mean = %.4f\n",
              p, mean(sigma_filt[p, ], na.rm = TRUE)))
}

# Clip to runner's allowed range
sigma_filt <- pmin(pmax(sigma_filt, 1e-6), 3.0)
ssm_sigma_per_plant <- rowMeans(sigma_filt, na.rm = TRUE)

# ============================================================================
# 3 . Direct: sigma_obs_pl[p] posterior median per plant.
# ============================================================================
direct_draws <- direct$draws
vars_avail   <- posterior::variables(direct_draws)
sigma_var <- if ("sigma_obs_pl[1]" %in% vars_avail) {
  paste0("sigma_obs_pl[", 1:P, "]")
} else if ("sigma_obs[1]" %in% vars_avail) {
  paste0("sigma_obs[", 1:P, "]")
} else {
  stop("Could not locate per-plant sigma parameter in Direct draws.")
}
direct_sigma_med <- apply(
  posterior::as_draws_matrix(
    posterior::subset_draws(direct_draws, variable = sigma_var)
  ),
  2, median
)

# ============================================================================
# 4 . Comparison table + verdict
# ============================================================================
comparison <- data.frame(
  plant            = plant_names,
  sigma_SSM_fixed  = round(ssm_sigma_per_plant, 4),
  sigma_Direct_med = round(direct_sigma_med, 4),
  ratio_Direct_SSM = round(direct_sigma_med / ssm_sigma_per_plant, 3),
  row.names = NULL
)

cat("\n============================================================\n")
cat("  sigma_obs partition: SSM (MARSS Stage 1) vs Direct (jointly estimated)\n")
cat("============================================================\n\n")
print(comparison, row.names = FALSE)

cat("\n[Aggregate]\n")
cat(sprintf("  Mean sigma_obs (SSM)    = %.4f\n", mean(ssm_sigma_per_plant)))
cat(sprintf("  Mean sigma_obs (Direct) = %.4f\n", mean(direct_sigma_med)))
cat(sprintf("  Mean ratio Direct/SSM   = %.3f\n",
            mean(direct_sigma_med / ssm_sigma_per_plant)))

out_csv <- file.path(PROJ_DIR, "Modelos AR con rho est. y rho=0",
                     "paper_pipeline", "results",
                     "sigma_obs_partition_SSM_vs_Direct.csv")
write.csv(comparison, out_csv, row.names = FALSE)
cat(sprintf("\n[CSV] wrote: %s\n", out_csv))

cat("\n[Verdict]\n")
n_lower <- sum(direct_sigma_med < ssm_sigma_per_plant)
if (n_lower >= 0.75 * P) {
  cat(sprintf("  Direct < SSM in %d/%d plants. Story HOLDS:\n", n_lower, P))
  cat("  Direct pulls sigma_obs down, pushing variance into rho.\n")
} else if (n_lower <= 0.25 * P) {
  cat(sprintf("  Direct > SSM in %d/%d plants. Story FAILS as stated;\n",
              P - n_lower, P))
  cat("  the shift in hat-rho is NOT explained by sigma_obs reassignment.\n")
} else {
  cat(sprintf("  Mixed: Direct < SSM in %d/%d, > SSM in %d/%d plants.\n",
              n_lower, P, P - n_lower, P))
  cat("  The naive story holds partially; use a cautious phrasing.\n")
}
cat("\n=== DONE ===\n")
