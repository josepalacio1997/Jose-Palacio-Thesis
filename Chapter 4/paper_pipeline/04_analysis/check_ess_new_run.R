# =============================================================================
# check_ess_new_run.R
# -----------------------------------------------------------------------------
# Compare ESS between the previous Direct fix run (stan_fits_fsv_v3, 3 chains x
# 8000 sampling = 24000 draws) and the new run (stan_fits_fsv_4ch9k, 4 chains x
# 9000 sampling = 36000 draws). The prediction was ess_min ~ 282 * 1.5 ~ 423.
#
# Run from the repo root:
#   Rscript "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/check_ess_new_run.R"
# =============================================================================

suppressPackageStartupMessages({
  library(posterior); library(dplyr)
})

PROJ_DIR <- "/Users/josepalacio/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW"

old_path <- file.path(PROJ_DIR, "outputs", "stan_fits_fsv_v3",    "rho_00.rds")
new_path <- file.path(PROJ_DIR, "outputs", "stan_fits_fsv_4ch9k", "rho_00.rds")

stopifnot(file.exists(old_path), file.exists(new_path))

cat("Loading fits (this takes ~30 seconds each)...\n")
old_fit <- readRDS(old_path)
new_fit <- readRDS(new_path)

ess_block <- function(fit, label) {
  draws <- fit$draws
  vars  <- posterior::variables(draws)
  # Block selectors: Rt[p,t], It[p,t], beta[p]
  rt_vars   <- grep("^Rt\\[",   vars, value = TRUE)
  it_vars   <- grep("^It\\[",   vars, value = TRUE)
  beta_vars <- grep("^beta\\[", vars, value = TRUE)

  block_stats <- function(var_names) {
    if (length(var_names) == 0) return(NULL)
    mat <- posterior::as_draws_matrix(
      posterior::subset_draws(draws, variable = var_names)
    )
    ess <- apply(mat, 2, posterior::ess_bulk)
    data.frame(
      ess_min = min(ess, na.rm = TRUE),
      ess_med = median(ess, na.rm = TRUE),
      ess_max = max(ess, na.rm = TRUE),
      n_params = length(ess)
    )
  }

  rbind(
    cbind(fit_label = label, block = "Rt",   block_stats(rt_vars)),
    cbind(fit_label = label, block = "It",   block_stats(it_vars)),
    cbind(fit_label = label, block = "beta", block_stats(beta_vars))
  )
}

old_stats <- ess_block(old_fit, "OLD (3ch x 8000)")
new_stats <- ess_block(new_fit, "NEW (4ch x 9000)")

cat("\n============================================================\n")
cat("  ESS comparison: Direct rho=0 — OLD vs NEW run\n")
cat("============================================================\n\n")

cmp <- rbind(old_stats, new_stats) |>
  arrange(block, fit_label) |>
  mutate(across(c(ess_min, ess_med, ess_max), \(x) round(x, 0)))
print(cmp, row.names = FALSE)

# Highlight the Rt min (the worst case from before, ess_min = 282)
old_rt_min <- old_stats$ess_min[old_stats$block == "Rt"]
new_rt_min <- new_stats$ess_min[new_stats$block == "Rt"]
cat("\n[Headline]\n")
cat(sprintf("  Rt ess_min:  OLD = %d  ->  NEW = %d  (predicted ~423)\n",
            old_rt_min, new_rt_min))
cat(sprintf("  Ratio NEW/OLD = %.2f x  (vs predicted 1.5 x)\n",
            new_rt_min / old_rt_min))

# Also check diagnostics
cat("\n[Diagnostics, NEW run]\n")
ds <- new_fit$diagnostic_summary
cat(sprintf("  iter_sampling (per chain): %d\n", new_fit$meta$iter_sampling))
cat(sprintf("  n_chains:                  %d\n", new_fit$meta$n_chains))
cat(sprintf("  num_divergent (sum):       %d\n", sum(ds$num_divergent)))
cat(sprintf("  num_max_treedepth (sum):   %d\n", sum(ds$num_max_treedepth)))
if (!is.null(ds$ebfmi))
  cat(sprintf("  EBFMI min:                 %.3f\n", min(ds$ebfmi)))

cat("\n=== DONE ===\n")
