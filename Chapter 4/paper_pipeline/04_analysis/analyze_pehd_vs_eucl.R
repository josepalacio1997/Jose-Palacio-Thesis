# =============================================================================
# analyze_pehd_vs_eucl.R
# -----------------------------------------------------------------------------
# Side-by-side comparison of the primary PEHD pipeline vs the Euclidean
# ablation. Loads the per-plant CRPS CSVs from both pipelines, produces
# a unified comparison table and figure, and reports differences in the
# posterior of rho between the two spatial weight matrices.
#
# Prerequisites:
#   - results_modelA_AR1/crps_by_plant_AR1.csv          (PEHD)
#   - results_modelA_AR1_eucl/crps_by_plant_AR1_eucl.csv (Euclidean)
#   - both pipelines saved a rho posterior summary table
#
# Outputs (saved to ../results/):
#   pehd_vs_eucl_crps.csv          : merged per-plant CRPS table
#   pehd_vs_eucl_rho.csv           : rho posterior summary, both pipelines
#   fig_pehd_vs_eucl_crps.pdf      : grouped bar chart (per-plant CRPS log10)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(ggplot2); library(scales)
})

PROJ_DIR <- "/Users/josepalacio/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW"
RESULTS_DIR <- file.path(PROJ_DIR, "Modelos AR con rho est. y rho=0",
                         "paper_pipeline", "results")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

PEHD_DIR <- file.path(PROJ_DIR, "Modelos AR con rho est. y rho=0",
                      "results_modelA_AR1")
EUCL_DIR <- file.path(PROJ_DIR, "Modelos AR con rho est. y rho=0",
                      "results_modelA_AR1_eucl")

# ============================================================================
# 1. Load per-plant CRPS from both pipelines
# ============================================================================
crps_pehd <- read_csv(file.path(PEHD_DIR, "crps_by_plant_AR1.csv"),
                      show_col_types = FALSE) |>
  mutate(pipeline = "PEHD")

crps_eucl <- read_csv(file.path(EUCL_DIR, "crps_by_plant_AR1_eucl.csv"),
                      show_col_types = FALSE) |>
  mutate(pipeline = "Euclidean")

# Keep only Direct variants (SSM not apples-to-apples vs y_raw)
crps_long <- bind_rows(crps_pehd, crps_eucl) |>
  filter(fit %in% c("A_direct_est", "A_direct_fix",
                    "A_joint_est",  "A_joint_fix")) |>  # accept old or new
  mutate(fit = recode(fit,
                      "A_joint_est" = "A_direct_est",
                      "A_joint_fix" = "A_direct_fix")) |>
  select(pipeline, fit, plant, n_obs, crps_mean, crps_log_mean)

# Wide form: PEHD vs Eucl, side by side
crps_wide <- crps_long |>
  pivot_wider(id_cols = c(fit, plant, n_obs),
              names_from = pipeline,
              values_from = c(crps_mean, crps_log_mean))

write_csv(crps_wide, file.path(RESULTS_DIR, "pehd_vs_eucl_crps.csv"))
cat("Per-plant comparison table:\n")
print(crps_wide, n = Inf)

# Aggregate
agg <- crps_long |>
  group_by(pipeline, fit) |>
  summarise(mean_crps     = mean(crps_mean,     na.rm = TRUE),
            mean_crps_log = mean(crps_log_mean, na.rm = TRUE),
            .groups = "drop")
cat("\nAggregate CRPS by pipeline and fit:\n")
print(agg)

# ============================================================================
# 2. Rho posterior comparison (PEHD vs Eucl, Direct and SSM)
# ============================================================================
rho_pehd <- read_csv(file.path(PEHD_DIR, "tableS_phi_posterior_AR1.csv"),
                     show_col_types = FALSE) |>
  mutate(pipeline = "PEHD")
rho_eucl <- read_csv(file.path(EUCL_DIR, "tableS_phi_posterior_AR1.csv"),
                     show_col_types = FALSE) |>
  mutate(pipeline = "Euclidean")

rho_tbl <- bind_rows(rho_pehd, rho_eucl) |>
  select(pipeline, fit, mean, sd, q5, q50, q95, rhat, ess_bulk)

write_csv(rho_tbl, file.path(RESULTS_DIR, "pehd_vs_eucl_rho.csv"))
cat("\nRho posterior summary, both pipelines:\n")
print(rho_tbl)

# ============================================================================
# 3. Figure: per-plant CRPS log10, PEHD vs Eucl, only Direct rho-est variant
# ============================================================================
plot_df <- crps_long |>
  filter(fit == "A_direct_est") |>
  mutate(plant = factor(plant, levels = paste0("WWTP", 1:15)))

p_crps <- ggplot(plot_df,
                 aes(x = plant, y = crps_log_mean, fill = pipeline)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.78) +
  scale_fill_manual(values = c("PEHD"      = "#2C7BB6",
                                "Euclidean" = "#D7191C")) +
  labs(x = NULL,
       y = expression("CRPS, log"[10]*" scale"),
       fill = "Spatial weight") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(RESULTS_DIR, "fig_pehd_vs_eucl_crps.pdf"),
       p_crps, width = 12, height = 5)

# ============================================================================
# 4. DONE
# ============================================================================
cat("\n=== DONE ===\n")
cat(sprintf("Wrote outputs to: %s\n", RESULTS_DIR))
cat("  pehd_vs_eucl_crps.csv\n")
cat("  pehd_vs_eucl_rho.csv\n")
cat("  fig_pehd_vs_eucl_crps.pdf\n")
