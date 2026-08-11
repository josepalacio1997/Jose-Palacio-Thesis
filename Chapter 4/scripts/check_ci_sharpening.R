# =============================================================
# Empirical check: does spatial coupling sharpen uncertainty
# at small plants via borrowing strength?
#
# Test: compare 95% CI widths of R_it and I_it between the
# rho_est and rho_fix variants, per plant, correlated with
# plant population.
# =============================================================

suppressPackageStartupMessages({
  library(posterior); library(dplyr); library(tidyr); library(ggplot2)
})

BASE     <- "~/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW/Modelos AR con rho est. y rho=0"
FITS_DIR <- file.path(BASE, "fits")
POP_PATH <- file.path(BASE, "../outputs/N_plantas.rds")
OUT_DIR  <- file.path(BASE, "results_modelA_AR1")

# ---- 1. Load fits ----
cat("Loading fits...\n")
fits <- list(
  Data_est    = readRDS(file.path(FITS_DIR, "A_direct_est.rds")),
  Data_fix    = readRDS(file.path(FITS_DIR, "A_direct_fix.rds")),
  Filter_est  = readRDS(file.path(FITS_DIR, "A_ssm_est.rds")),
  Filter_fix  = readRDS(file.path(FITS_DIR, "A_ssm_fix.rds"))
)

N_plantas <- readRDS(POP_PATH)
plant_names <- paste0("WWTP", seq_along(N_plantas))

# ---- 2. Helper: get 95% CI width per (plant, week) ----
ci_width_by_plant <- function(save_obj, param) {
  d <- as_draws_matrix(save_obj$draws)
  vars <- grep(paste0("^", param, "\\["), colnames(d), value = TRUE)
  q <- apply(d[, vars, drop = FALSE], 2,
             quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  width <- q[2, ] - q[1, ]

  parsed <- do.call(rbind,
    lapply(strsplit(sub(".*\\[(.*)\\]", "\\1", vars), ","),
           function(x) as.integer(x)))
  out <- data.frame(
    plant_idx = parsed[, 1],
    week      = parsed[, 2],
    width     = width
  )
  # Average over weeks per plant
  agg <- out |>
    group_by(plant_idx) |>
    summarise(mean_width = mean(width, na.rm = TRUE),
              median_width = median(width, na.rm = TRUE),
              .groups = "drop") |>
    arrange(plant_idx)
  agg$plant <- plant_names[agg$plant_idx]
  agg
}

# ---- 3. Compute CI widths for Rt and It, all fits ----
cat("Computing CI widths (this takes ~1 min per fit)...\n")
widths <- list()
for (nm in names(fits)) {
  cat(sprintf("  %s Rt ... ", nm))
  w_rt <- ci_width_by_plant(fits[[nm]], "Rt")
  cat("It ... ")
  w_it <- ci_width_by_plant(fits[[nm]], "It")
  cat("done\n")
  widths[[nm]] <- list(Rt = w_rt, It = w_it)
}

# ---- 4. Build comparison table ----
build_compare <- function(param) {
  d_est <- widths$Data_est[[param]]   |> select(plant, mean_width) |> rename(width_data_est   = mean_width)
  d_fix <- widths$Data_fix[[param]]   |> select(plant, mean_width) |> rename(width_data_fix   = mean_width)
  f_est <- widths$Filter_est[[param]] |> select(plant, mean_width) |> rename(width_filter_est = mean_width)
  f_fix <- widths$Filter_fix[[param]] |> select(plant, mean_width) |> rename(width_filter_fix = mean_width)

  cmp <- Reduce(function(a, b) full_join(a, b, by = "plant"),
                list(d_est, d_fix, f_est, f_fix)) |>
    mutate(
      pop_K            = N_plantas[plant] / 1000,
      ratio_data       = width_data_est   / width_data_fix,
      ratio_filter     = width_filter_est / width_filter_fix,
      param            = param
    ) |>
    arrange(pop_K)
  cmp
}

rt_cmp <- build_compare("Rt")
it_cmp <- build_compare("It")

cat("\n===== Rt: 95% CI width comparison =====\n")
print(rt_cmp |> select(plant, pop_K, width_data_est, width_data_fix,
                        ratio_data, width_filter_est, width_filter_fix,
                        ratio_filter),
      row.names = FALSE, digits = 3)

cat("\n===== It: 95% CI width comparison =====\n")
print(it_cmp |> select(plant, pop_K, width_data_est, width_data_fix,
                        ratio_data, width_filter_est, width_filter_fix,
                        ratio_filter),
      row.names = FALSE, digits = 3)

# ---- 5. Correlations ----
cat("\n===== Correlations: CI-width RATIO vs population =====\n")
cat("Hypothesis: negative correlation (small plants shrink more)\n\n")
cat(sprintf("Rt Data   ratio vs pop: r = %.3f (Pearson), rho = %.3f (Spearman)\n",
            cor(rt_cmp$ratio_data, rt_cmp$pop_K, method = "pearson", use = "complete.obs"),
            cor(rt_cmp$ratio_data, rt_cmp$pop_K, method = "spearman", use = "complete.obs")))
cat(sprintf("Rt Filter ratio vs pop: r = %.3f (Pearson), rho = %.3f (Spearman)\n",
            cor(rt_cmp$ratio_filter, rt_cmp$pop_K, method = "pearson", use = "complete.obs"),
            cor(rt_cmp$ratio_filter, rt_cmp$pop_K, method = "spearman", use = "complete.obs")))
cat(sprintf("It Data   ratio vs pop: r = %.3f (Pearson), rho = %.3f (Spearman)\n",
            cor(it_cmp$ratio_data, it_cmp$pop_K, method = "pearson", use = "complete.obs"),
            cor(it_cmp$ratio_data, it_cmp$pop_K, method = "spearman", use = "complete.obs")))
cat(sprintf("It Filter ratio vs pop: r = %.3f (Pearson), rho = %.3f (Spearman)\n",
            cor(it_cmp$ratio_filter, it_cmp$pop_K, method = "pearson", use = "complete.obs"),
            cor(it_cmp$ratio_filter, it_cmp$pop_K, method = "spearman", use = "complete.obs")))

# ---- 6. Save CSVs and a summary plot ----
write.csv(rt_cmp, file.path(OUT_DIR, "ci_width_Rt_by_plant.csv"), row.names = FALSE)
write.csv(it_cmp, file.path(OUT_DIR, "ci_width_It_by_plant.csv"), row.names = FALSE)
cat("\n  -> ci_width_Rt_by_plant.csv, ci_width_It_by_plant.csv\n")

# Scatter plot: ratio vs population
plot_df <- bind_rows(
  rt_cmp |> select(plant, pop_K, ratio_data, ratio_filter) |>
    pivot_longer(c(ratio_data, ratio_filter), names_to = "variant", values_to = "ratio") |>
    mutate(param = "R[it]"),
  it_cmp |> select(plant, pop_K, ratio_data, ratio_filter) |>
    pivot_longer(c(ratio_data, ratio_filter), names_to = "variant", values_to = "ratio") |>
    mutate(param = "I[it]")
) |>
  mutate(variant = recode(variant, "ratio_data" = "Data", "ratio_filter" = "Filter"))

p <- ggplot(plot_df, aes(x = pop_K, y = ratio, color = variant, label = plant)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
  facet_wrap(~ param, labeller = label_parsed) +
  scale_x_log10() +
  scale_color_manual(values = c("Data" = "steelblue", "Filter" = "darkorange")) +
  labs(x = "Population (K, log scale)",
       y = "CI width ratio (rho_est / rho_fix)",
       title = "Does spatial coupling narrow CIs at small plants?",
       subtitle = "Ratio < 1: coupling sharpens; > 1: coupling widens") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(OUT_DIR, "ci_width_ratio_vs_pop.pdf"),
       p, width = 10, height = 5)
cat("  -> ci_width_ratio_vs_pop.pdf\n")

cat("\nDone.\n")
