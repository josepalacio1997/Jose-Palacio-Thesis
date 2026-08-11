# =============================================================================
# diagnose_chain_disagreement.R   (v2 — correct ESS via fit$summary)
# -----------------------------------------------------------------------------
# Recompute per-parameter ESS comparing OLD and NEW Direct (rho=0) runs.
# Uses the fit$summary table (cmdstanr's per-variable ess_bulk and rhat),
# NOT manual posterior::ess_bulk() calls on as_draws_matrix output (which
# silently aggregate across variables and give wrong scalar values).
#
# Also: if ESS truly collapsed, identify the outlier chain by computing
# per-chain posterior means for the worst-ess Rt parameters.
#
# Run from the repo root:
#   Rscript "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/diagnose_chain_disagreement.R"
# =============================================================================

suppressPackageStartupMessages({
  library(posterior); library(dplyr); library(stringr)
})

PROJ_DIR <- "/Users/josepalacio/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW"
old_path <- file.path(PROJ_DIR, "outputs", "stan_fits_fsv_v3",    "rho_00.rds")
new_path <- file.path(PROJ_DIR, "outputs", "stan_fits_fsv_4ch9k", "rho_00.rds")
stopifnot(file.exists(old_path), file.exists(new_path))

old_fit <- readRDS(old_path)
new_fit <- readRDS(new_path)

# fit$summary has per-variable ess_bulk and rhat from the cmdstanr fit.
# Each row is one parameter.
block_stats <- function(fit, label) {
  s <- fit$summary
  # Tag block
  s$block <- dplyr::case_when(
    grepl("^Rt\\[",   s$variable) ~ "Rt",
    grepl("^It\\[",   s$variable) ~ "It",
    grepl("^beta\\[", s$variable) ~ "beta",
    grepl("^sigma_obs_pl\\[", s$variable) ~ "sigma_obs_pl",
    TRUE                          ~ NA_character_
  )
  s |>
    dplyr::filter(!is.na(block)) |>
    dplyr::group_by(block) |>
    dplyr::summarise(
      n_params = n(),
      ess_min  = min(ess_bulk, na.rm = TRUE),
      ess_med  = median(ess_bulk, na.rm = TRUE),
      ess_max  = max(ess_bulk, na.rm = TRUE),
      rhat_max = max(rhat,    na.rm = TRUE),
      rhat_med = median(rhat, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(fit_label = label, .before = 1)
}

old_stats <- block_stats(old_fit, "OLD (3ch x 8000)")
new_stats <- block_stats(new_fit, "NEW (4ch x 9000)")

cat("============================================================\n")
cat("  ESS/Rhat per block — OLD vs NEW (using fit$summary, correct)\n")
cat("============================================================\n\n")
print(dplyr::bind_rows(old_stats, new_stats) |> dplyr::arrange(block, fit_label),
      row.names = FALSE)

# ---- Headline contrast (worst-case Rt ess_min) ------------------------------
old_rt <- old_stats$ess_min[old_stats$block == "Rt"]
new_rt <- new_stats$ess_min[new_stats$block == "Rt"]
cat(sprintf("\n[Headline] Rt ess_min:  OLD = %.0f -> NEW = %.0f  (predicted ~423)\n",
            old_rt, new_rt))
cat(sprintf("           ratio NEW/OLD = %.2fx\n", new_rt / old_rt))

# ============================================================================
# If NEW Rt ess_min < OLD: investigate which chain disagrees.
# ============================================================================
if (new_rt < old_rt) {
  cat("\n*** NEW ess_min is WORSE than OLD. Investigating chain disagreement ***\n")
  # Top-5 worst Rt parameters in the new fit
  s_new <- new_fit$summary
  rt_new <- s_new[grepl("^Rt\\[", s_new$variable), ]
  worst5 <- rt_new[order(rt_new$ess_bulk)[1:5], c("variable", "ess_bulk", "rhat")]
  cat("\n[Top 5 worst Rt parameters in NEW]\n")
  print(worst5, row.names = FALSE)

  # Per-chain stats for the worst variable
  worst_var <- worst5$variable[1]
  cat(sprintf("\n[Per-chain stats for worst: %s]\n", worst_var))
  draws_arr <- new_fit$draws
  mat <- posterior::as_draws_matrix(
    posterior::subset_draws(draws_arr, variable = worst_var)
  )
  n_iter  <- dim(draws_arr)[1]
  n_chain <- dim(draws_arr)[2]
  cid <- rep(1:n_chain, each = n_iter)
  per_chain <- data.frame(
    chain = 1:n_chain,
    mean  = tapply(as.numeric(mat), cid, mean),
    sd    = tapply(as.numeric(mat), cid, sd),
    q05   = tapply(as.numeric(mat), cid, quantile, probs = 0.05),
    q95   = tapply(as.numeric(mat), cid, quantile, probs = 0.95)
  )
  print(per_chain, row.names = FALSE)
  cat(sprintf("\n  Inter-chain spread of means: max - min = %.4g\n",
              max(per_chain$mean) - min(per_chain$mean)))
  cat(sprintf("  Median within-chain SD:       %.4g\n",
              median(per_chain$sd)))
  cat(sprintf("  Ratio inter/within:            %.2g  (>1 = chains disagree more than they vary internally)\n",
              (max(per_chain$mean) - min(per_chain$mean)) / median(per_chain$sd)))
} else {
  cat("\n[OK] NEW Rt ess_min >= OLD; no chain disagreement to investigate.\n")
}

cat("\n=== DONE ===\n")
