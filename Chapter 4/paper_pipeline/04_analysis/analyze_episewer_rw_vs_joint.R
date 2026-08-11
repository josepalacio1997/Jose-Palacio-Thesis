# =============================================================================
# analyze_episewer_rw_vs_joint.R
# -----------------------------------------------------------------------------
# Compares EpiSewer DAILY+RW (single-plant baseline) against the proposed
# Direct first-order model on CRPS, both on linear (B gc/day) and log10 scales.
#
# EpiSewer RW pipeline:
#   - input concentration = y_raw (B gc/day), flow = 1 (placeholder)
#   - predicted_concentration is DIRECTLY in B gc/day (no conversion needed)
#   - input_data carries real_date + pseudo_date for downstream mapping
#
# Output:
#   crps_episewer_rw_vs_direct.csv
#   fig_episewer_rw_vs_direct.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(posterior); library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(scoringRules); library(readr)
})

PROJ_DIR    <- "/Users/josepalacio/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW"
EP_FITS_DIR <- file.path(PROJ_DIR, "Modelos AR con rho est. y rho=0",
                          "episewer_rw_fits")
OUT_DIR     <- file.path(PROJ_DIR, "Modelos AR con rho est. y rho=0",
                          "results_modelA_AR1")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

P <- 15
plant_names <- paste0("WWTP", 1:P)

# ============================================================================
# 1 . Para cada planta, computar CRPS de EpiSewer RW
#     predicted_concentration ya esta en B gc/day -> sin conversion
# ============================================================================
cat("============================================================\n")
cat("  Computing EpiSewer RW CRPS per plant\n")
cat("============================================================\n")
episewer_crps <- list()
for (plant_idx in 1:15) {
  plant_name <- sprintf("WWTP%d", plant_idx)
  rds_path <- file.path(EP_FITS_DIR, sprintf("plant_%02d.rds", plant_idx))
  if (!file.exists(rds_path)) {
    cat(sprintf("  %s: MISSING\n", plant_name)); next
  }
  cat(sprintf("  %s ... ", plant_name))
  t0 <- Sys.time()
  ep <- tryCatch(readRDS(rds_path), error = function(e) NULL)
  if (is.null(ep)) {cat("READ ERROR\n"); next}

  # Guard: si EpiSewer descarto el fit -> skip
  fit <- ep$result$fit
  if (is.null(fit) || !inherits(fit, "CmdStanMCMC")) {
    cat(sprintf("SKIP (no fit, size=%.0f KB)\n", file.size(rds_path)/1024))
    next
  }

  # y_obs = measurement (B gc/day) del input_data
  y_obs       <- ep$input_data$measurement
  real_dates  <- as.Date(ep$input_data$real_date)
  N_obs <- length(y_obs)

  # Predictive draws: predicted_concentration
  # Las primeras N_obs columnas corresponden a las obs (resto = forecast horizon)
  pc_draws <- fit$draws("predicted_concentration", format = "draws_matrix")
  N_pred <- ncol(pc_draws)
  if (N_pred < N_obs) {
    cat(sprintf("FAIL (n_pred=%d < n_obs=%d)\n", N_pred, N_obs)); next
  }

  # Tomamos solo las primeras N_obs columnas (en sample)
  pc_obs <- pc_draws[, seq_len(N_obs), drop = FALSE]

  # CRPS contra y_obs (ya estan en mismas unidades B gc/day)
  crps_orig <- rep(NA_real_, N_obs)
  crps_log  <- rep(NA_real_, N_obs)
  for (i in seq_len(N_obs)) {
    yo <- y_obs[i]
    yp <- pc_obs[, i]
    yp_pos <- yp[is.finite(yp) & yp > 0]
    if (is.na(yo) || yo <= 0 || length(yp_pos) < 20) next
    crps_orig[i] <- scoringRules::crps_sample(y = yo, dat = yp_pos)
    crps_log[i]  <- scoringRules::crps_sample(y = log10(yo),
                                                dat = log10(yp_pos))
  }

  episewer_crps[[plant_name]] <- data.frame(
    plant                   = plant_name,
    n_obs                   = sum(!is.na(crps_orig)),
    EpiSewer_crps_mean      = mean(crps_orig, na.rm = TRUE),
    EpiSewer_crps_median    = median(crps_orig, na.rm = TRUE),
    EpiSewer_crps_log_mean  = mean(crps_log,  na.rm = TRUE),
    EpiSewer_crps_log_median= median(crps_log, na.rm = TRUE)
  )
  cat(sprintf("CRPS_mean=%.2f, CRPS_log=%.4f (%.1f sec)\n",
              episewer_crps[[plant_name]]$EpiSewer_crps_mean,
              episewer_crps[[plant_name]]$EpiSewer_crps_log_mean,
              as.numeric(difftime(Sys.time(), t0, units="secs"))))
}

ep_summary <- do.call(rbind, episewer_crps); rownames(ep_summary) <- NULL
cat(sprintf("\nEpiSewer RW summary (%d plants):\n", nrow(ep_summary)))
print(ep_summary, row.names = FALSE, digits = 4)

# ============================================================================
# 2 . Cargar Joint CRPS by plant (de crps_by_plant_AR1.csv si existe,
#     o el global crps_by_plant.csv)
# ============================================================================
crps_csv_AR1 <- file.path(OUT_DIR, "crps_by_plant_AR1.csv")
crps_csv_old <- file.path(PROJ_DIR, "Modelos AR con rho est. y rho=0",
                          "results_modelA", "crps_by_plant.csv")
crps_csv <- if (file.exists(crps_csv_AR1)) crps_csv_AR1 else crps_csv_old
cat(sprintf("\nLeyendo Joint CRPS de: %s\n", basename(crps_csv)))

# NOTE: SSM omitted — SSM is trained against y_filt (MARSS Stage-1 output),
# while Direct and EpiSewer are trained against y_raw. Scoring SSM against
# y_raw is not apples-to-apples. Only Direct vs EpiSewer is a fair comparison.
direct_by_plant <- read.csv(crps_csv) |>
  filter(fit %in% c("A_direct_est", "A_direct_fix",
                    "A_joint_est",  "A_joint_fix")) |>  # accept old or new
  mutate(fit = recode(fit,
                      "A_joint_est" = "A_direct_est",
                      "A_joint_fix" = "A_direct_fix")) |>
  select(fit, plant, crps_mean, crps_log_mean) |>
  pivot_wider(names_from = fit,
              values_from = c(crps_mean, crps_log_mean),
              names_glue = "{fit}_{.value}")

# ============================================================================
# 3 . Merge + ratios
# ============================================================================
compare <- ep_summary |>
  left_join(direct_by_plant, by = "plant") |>
  mutate(
    ratio_orig_vs_direct_est = A_direct_est_crps_mean      / EpiSewer_crps_mean,
    ratio_log_vs_direct_est  = A_direct_est_crps_log_mean  / EpiSewer_crps_log_mean,
    ratio_orig_vs_direct_fix = A_direct_fix_crps_mean      / EpiSewer_crps_mean,
    ratio_log_vs_direct_fix  = A_direct_fix_crps_log_mean  / EpiSewer_crps_log_mean
  )

cat("\n--- compare table (Direct only — SSM omitted, uses y_filt) ---\n")
print(compare |> select(plant, n_obs, EpiSewer_crps_log_mean,
                         A_direct_est_crps_log_mean, A_direct_fix_crps_log_mean,
                         ratio_log_vs_direct_est, ratio_log_vs_direct_fix),
      row.names = FALSE, digits = 3)

write.csv(compare, file.path(OUT_DIR, "crps_episewer_rw_vs_direct.csv"),
          row.names = FALSE)
cat(sprintf("\n  -> %s\n", file.path(OUT_DIR, "crps_episewer_rw_vs_direct.csv")))

# ============================================================================
# 4 . Agregado global
# ============================================================================
cat("\n=== AGGREGATE across all plants with data ===\n")
agg <- compare |>
  summarise(
    n_plants             = n(),
    EpiSewer_mean        = mean(EpiSewer_crps_mean,         na.rm=TRUE),
    EpiSewer_log_mean    = mean(EpiSewer_crps_log_mean,     na.rm=TRUE),
    Direct_est_mean      = mean(A_direct_est_crps_mean,     na.rm=TRUE),
    Direct_est_log_mean  = mean(A_direct_est_crps_log_mean, na.rm=TRUE),
    Direct_fix_mean      = mean(A_direct_fix_crps_mean,     na.rm=TRUE),
    Direct_fix_log_mean  = mean(A_direct_fix_crps_log_mean, na.rm=TRUE)
  )
print(agg, digits = 4)

cat(sprintf("\n  Direct_est vs EpiSewer: %d/%d plants where Direct_est is better (ratio_log < 1)\n",
            sum(compare$ratio_log_vs_direct_est < 1, na.rm=TRUE), nrow(compare)))
cat(sprintf("  Direct_fix vs EpiSewer: %d/%d plants where Direct_fix is better\n",
            sum(compare$ratio_log_vs_direct_fix < 1, na.rm=TRUE), nrow(compare)))

# ============================================================================
# 5 . Figure: per-plant log10 CRPS — Direct first-order vs EpiSewer baseline
# ============================================================================
long <- compare |>
  select(plant,
         EpiSewer        = EpiSewer_crps_log_mean,
         `Direct; rho-est` = A_direct_est_crps_log_mean,
         `Direct; rho=0`   = A_direct_fix_crps_log_mean) |>
  pivot_longer(-plant, names_to = "Model", values_to = "CRPS_log10") |>
  mutate(plant = factor(plant, levels = paste0("WWTP", 1:15)),
         Model = factor(Model,
                        levels = c("Direct; rho-est",
                                   "Direct; rho=0",
                                   "EpiSewer")))

# Labels typeset with plotmath: "Direct; rho-hat", "Direct; rho == 0", "EpiSewer"
legend_labels <- c(
  "Direct; rho-est" = expression("Direct; " * hat(rho)),
  "Direct; rho=0"   = expression("Direct; " * rho == 0),
  "EpiSewer"          = expression("EpiSewer")
)

p <- ggplot(long, aes(x = plant, y = CRPS_log10, fill = Model)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.78) +
  scale_fill_manual(values = c("Direct; rho-est" = "#2C7BB6",
                                "Direct; rho=0"   = "#9ECAE1",
                                "EpiSewer"          = "#5E5E5E"),
                    labels = legend_labels) +
  labs(x = NULL, y = expression("CRPS, log"[10]*" scale"), fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(OUT_DIR, "fig_episewer_rw_vs_direct.pdf"),
       p, width = 14, height = 5.5)
cat(sprintf("\n  -> %s\n", file.path(OUT_DIR, "fig_episewer_rw_vs_direct.pdf")))

cat("\n=== DONE ===\n")
