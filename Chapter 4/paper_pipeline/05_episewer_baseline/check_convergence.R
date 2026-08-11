# =============================================================================
# check_convergence_subset.R
# -----------------------------------------------------------------------------
# Convergencia para subset de plantas (por defecto las 3 que tenemos).
# Reporta: size, rhat_R, rhat_pred, ESS_R, ESS_min, divergences, treedepth,
# EBFMI, y veredicto (OK/WARN/FAIL).
# =============================================================================

suppressPackageStartupMessages({
  library(posterior); library(dplyr)
})

PROJ_DIR <- "/Users/josepalacio/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW"
EP_DIR   <- file.path(PROJ_DIR, "Modelos AR con rho est. y rho=0",
                      "episewer_rw_fits")

# Plantas a verificar
# Auto-detecta todas las plantas disponibles
files <- list.files(EP_DIR, pattern = "plant_\\d+\\.rds")
PLANT_IDX <- sort(as.integer(sub("plant_(\\d+)\\.rds", "\\1", files)))

cat("============================================================\n")
cat("  EpiSewer RW convergence subset\n")
cat("============================================================\n\n")

results <- list()
for (n in PLANT_IDX) {
  rds <- file.path(EP_DIR, sprintf("plant_%02d.rds", n))
  if (!file.exists(rds)) {
    cat(sprintf("plant_%02d: MISSING\n", n)); next
  }
  cat(sprintf("plant_%02d ... ", n))
  t0 <- Sys.time()
  ep <- readRDS(rds)
  fit <- ep$result$fit
  if (is.null(fit) || !inherits(fit, "CmdStanMCMC")) {
    cat("NO FIT\n"); next
  }

  size_MB <- file.size(rds) / 1e6
  iters <- if (!is.null(ep$meta$iter_sampling)) ep$meta$iter_sampling else NA
  diag <- fit$diagnostic_summary(quiet = TRUE)
  sm <- fit$summary()

  rhat_max  <- max(sm$rhat, na.rm = TRUE)
  rhat_R    <- max(sm$rhat[grepl("^R\\[", sm$variable)], na.rm = TRUE)
  rhat_pred <- max(sm$rhat[grepl("^predicted_concentration\\[", sm$variable)], na.rm = TRUE)
  ess_min   <- min(sm$ess_bulk, na.rm = TRUE)
  ess_R     <- min(sm$ess_bulk[grepl("^R\\[", sm$variable)], na.rm = TRUE)
  n_div     <- sum(diag$num_divergent)
  n_tree    <- sum(diag$num_max_treedepth)
  ebfmi_min <- min(diag$ebfmi, na.rm = TRUE)

  # Veredicto
  total_iters <- 4 * iters   # 4 chains × iter_sampling
  pct_div <- n_div / total_iters
  verdict <- "OK"
  if (rhat_R > 1.05)         verdict <- "FAIL"
  else if (rhat_pred > 1.05) verdict <- "FAIL"
  else if (ess_R < 400)      verdict <- "WARN (low ESS_R)"
  else if (pct_div > 0.01)   verdict <- "WARN (>1% div)"

  results[[as.character(n)]] <- data.frame(
    plant = sprintf("WWTP%d", n),
    size_MB = round(size_MB, 1),
    iters = iters,
    rhat_max = round(rhat_max, 4),
    rhat_R = round(rhat_R, 4),
    rhat_pred = round(rhat_pred, 4),
    ESS_min = round(ess_min),
    ESS_R = round(ess_R),
    div = n_div,
    treedepth = n_tree,
    ebfmi_min = round(ebfmi_min, 3),
    verdict = verdict
  )

  cat(sprintf("rhat_R=%.4f, rhat_pred=%.4f, ESS_R=%.0f, div=%d  [%s] (%.0f sec)\n",
              rhat_R, rhat_pred, ess_R, n_div, verdict,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

cat("\n============================================================\n")
cat("  Tabla completa\n")
cat("============================================================\n")
df <- do.call(rbind, results); rownames(df) <- NULL
print(df, row.names = FALSE)

# Veredicto final
cat("\n=== Resumen ===\n")
print(table(df$verdict))
