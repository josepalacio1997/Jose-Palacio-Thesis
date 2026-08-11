# =============================================================================
# regenerate_filter_estimates.R
# -----------------------------------------------------------------------------
# Regenerates the filter_estimates_15.pdf plot WITH gray bars marking weeks
# where the raw y_it is missing (which the original QMD version was missing).
#
# Uses y_raw, y_filt, sigma_filt already stored inside the A_ssm_est fit.
# Output: filter_estimates_15_v2.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

BASE     <- "~/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW/Modelos AR con rho est. y rho=0"
FIT_PATH <- file.path(BASE, "fits/A_ssm_est.rds")
POP_CACHE <- file.path(BASE, "../outputs/N_plantas.rds")
OUT_PATH  <- file.path(BASE, "..", "filter_estimates_15_v2.pdf")

# ---- 1. Load fit (contains y_raw + y_filt + sigma_filt) ----
cat("Loading fit...\n")
obj <- readRDS(FIT_PATH)
y_raw      <- obj$y_raw        # [P x T], NA where missing
y_filt     <- obj$y_filt       # [P x T], MARSS-smoothed
sigma_filt <- obj$sigma_filt   # [P x T]
P <- nrow(y_raw); T_wks <- ncol(y_raw)
cat(sprintf("  %d plants x %d weeks; %d missing y_raw (%.1f%%)\n",
            P, T_wks, sum(is.na(y_raw)), 100*mean(is.na(y_raw))))

# ---- 2. Plant labels with populations ----
plant_names <- paste0("WWTP", seq_len(P))
if (file.exists(POP_CACHE)) {
  N_plantas <- readRDS(POP_CACHE)
  plant_labels <- setNames(
    ifelse(is.na(N_plantas),
           plant_names,
           sprintf("%s (%dK)", plant_names, round(N_plantas / 1000))),
    plant_names
  )
} else {
  plant_labels <- setNames(plant_names, plant_names)
}

# ---- 3. Build long-format data ----
plot_df <- do.call(rbind, lapply(seq_len(P), function(p) {
  psi   <- y_filt[p, ]
  sig   <- sigma_filt[p, ]
  mu_ln <- log(pmax(psi, 1e-300)) - 0.5 * sig^2
  data.frame(
    week     = seq_len(T_wks),
    plant    = plant_names[p],
    observed = as.numeric(y_raw[p, ]),
    filtered = psi,
    lo       = exp(mu_ln + qnorm(0.025) * sig),
    hi       = exp(mu_ln + qnorm(0.975) * sig)
  )
}))
plot_df$plant <- factor(plot_df$plant, levels = plant_names,
                        labels = plant_labels[plant_names])

# ---- 4. Missing-week mask (for gray bars) ----
miss_df <- do.call(rbind, lapply(seq_len(P), function(p) {
  data.frame(
    week  = which(is.na(y_raw[p, ])),
    plant = plant_names[p]
  )
}))
miss_df$plant <- factor(miss_df$plant, levels = plant_names,
                        labels = plant_labels[plant_names])
cat(sprintf("  gray bars will mark %d missing weeks across all plants\n",
            nrow(miss_df)))

# ---- 5. Plot ----
p_filter <- ggplot(plot_df, aes(x = week)) +
  # Gray bars for missing raw y_it weeks
  geom_rect(data = miss_df,
            aes(xmin = week - 0.5, xmax = week + 0.5,
                ymin = -Inf, ymax = Inf),
            fill = "grey85", alpha = 0.5, inherit.aes = FALSE) +
  # Filter CI band
  geom_ribbon(aes(ymin = lo, ymax = hi),
              fill = "steelblue", alpha = 0.25, na.rm = TRUE) +
  # Filter mean
  geom_line(aes(y = filtered), colour = "steelblue4",
            linewidth = 0.6, na.rm = TRUE) +
  # Observed points (skips NAs)
  geom_point(aes(y = observed), colour = "black",
             size = 0.7, na.rm = TRUE) +
  facet_wrap(~ plant, ncol = 3, nrow = 5, scales = "free_y") +
  labs(x = "Week", y = "Viral load (B gc/day)") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92"),
        panel.grid.minor = element_blank(),
        strip.text       = element_text(face = "bold"))

# ---- 6. Save ----
ggsave(OUT_PATH, p_filter, width = 12, height = 8)
cat("\n  -> ", normalizePath(OUT_PATH, mustWork = FALSE), "\n", sep = "")
cat("\nDone.\n")
