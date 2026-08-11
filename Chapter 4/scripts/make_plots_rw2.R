# =============================================================
# Panel plots for the second-order random walk (RW2) fits
#   - Second-order random walk on lambda: x_t = 2*x_{t-1} - x_{t-2} + w_t
#     (not to be confused with AR(2), where coefficients are estimated)
#   - Filter (ssm) variants overlay y_filt (Stage-1 smoothed)
#   - Data (direct) variants overlay y_raw (observed viral load)
#   - Gray bars mark weeks where raw y_it is missing
#   - Rt plots include a dashed horizontal line at R = 1
#
# NOTE: only rho-estimated fits exist for RW2 (ssm_fix and direct_fix
# were not generated). Script produces 4 PDFs:
#   Rit_filter_rho_est_rw2.pdf,  Iit_filter_rho_est_rw2.pdf
#   Rit_data_rho_est_rw2.pdf,    Iit_data_rho_est_rw2.pdf
#
# Fit filenames still use the legacy `_ar2` suffix internally
# (A_ssm_est_ar2.rds, A_direct_est_ar2.rds) because that's what
# the earlier fitting scripts wrote to disk. Output PDFs use `rw2`
# for the audience-facing correct terminology.
# =============================================================

suppressPackageStartupMessages({
  library(posterior)
  library(ggplot2); library(dplyr); library(tidyr)
  library(ggh4x)
})

# Colors: Rt in orange, It in blue
COL_RT <- "darkorange"
COL_IT <- "steelblue"

# ---- PATHS ----
BASE <- "~/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW/Modelos AR con rho est. y rho=0"
# NOTE: fit files still named _ar2 on disk; content is RW2
FIT_SSM_EST_RW2    <- file.path(BASE, "fits/A_ssm_est_ar2.rds")
FIT_DIRECT_EST_RW2 <- file.path(BASE, "fits/A_direct_est_ar2.rds")
OUT_DIR            <- file.path(BASE, "results_modelA_RW2")
POP_CACHE          <- file.path(BASE, "../outputs/N_plantas.rds")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load fits (contain draws + y_raw + y_filt) ----
obj_ssm    <- readRDS(FIT_SSM_EST_RW2)
obj_direct <- readRDS(FIT_DIRECT_EST_RW2)
cat("Loaded RW2 fits.\n")

y_raw  <- obj_ssm$y_raw
y_filt <- obj_ssm$y_filt

m <- nrow(y_raw); T_wks <- ncol(y_raw)
miss_mask <- is.na(y_raw)
cat(sprintf("Data: %d plants x %d weeks; %d missing (%.1f%%)\n",
            m, T_wks, sum(miss_mask), 100*mean(miss_mask)))

# ---- Week -> Date mapping (uses actual dates from fit; falls back to guess) --
week_dates <- suppressWarnings(as.Date(colnames(y_raw)))
if (all(is.na(week_dates))) {
  warning("colnames(y_raw) not parseable as dates; using START_DATE = 2023-05-01")
  week_dates <- seq(from = as.Date("2023-05-01"),
                    by = "1 week", length.out = T_wks)
}
cat(sprintf("Date range: %s -- %s\n",
            format(min(week_dates), "%Y-%m-%d"),
            format(max(week_dates), "%Y-%m-%d")))

# ---- 1b. Plant labels with populations ----
plant_names <- paste0("WWTP", seq_len(m))
if (file.exists(POP_CACHE)) {
  N_plantas <- readRDS(POP_CACHE)
  plant_labels <- setNames(
    ifelse(is.na(N_plantas),
           plant_names,
           sprintf("%s (%dK)", plant_names, round(N_plantas / 1000))),
    plant_names
  )
} else {
  warning("N_plantas.rds not found -- falling back to bare WWTP labels")
  plant_labels <- setNames(plant_names, plant_names)
}

# ---- 2. Extract posterior summaries for Rt and It ----
extract_post <- function(save_obj, param) {
  d <- as_draws_matrix(save_obj$draws)
  vars <- grep(paste0("^", param, "\\["), colnames(d), value = TRUE)
  q <- apply(d[, vars, drop = FALSE], 2,
             quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  out <- data.frame(var = vars,
                    lo = q[1, ], md = q[2, ], hi = q[3, ])
  parsed <- do.call(rbind,
    lapply(strsplit(sub(".*\\[(.*)\\]", "\\1", vars), ","),
           function(x) as.integer(x)))
  out$plant_idx <- parsed[, 1]
  out$week      <- parsed[, 2]
  out$plant     <- plant_names[out$plant_idx]
  out
}

cat("Extracting Rt and It posteriors...\n")
Rt_ssm    <- extract_post(obj_ssm,    "Rt")
Rt_direct <- extract_post(obj_direct, "Rt")
It_ssm    <- extract_post(obj_ssm,    "It")
It_direct <- extract_post(obj_direct, "It")
cat("  done.\n")

# ---- 3. Plot function ----
make_plot <- function(post_df, overlay_mat, plants_subset,
                      y_lab, right_lab, ribbon_color = "steelblue",
                      add_hline_1 = FALSE) {
  if (is.null(rownames(overlay_mat))) rownames(overlay_mat) <- plant_names
  if (is.null(rownames(miss_mask)))   rownames(miss_mask)   <- plant_names

  df <- post_df[post_df$plant %in% plants_subset, ]
  df$plant <- factor(df$plant, levels = plants_subset,
                     labels = plant_labels[plants_subset])
  df$date  <- week_dates[df$week]

  ov <- data.frame(
    plant = rep(plants_subset, each = T_wks),
    week  = rep(seq_len(T_wks), times = length(plants_subset)),
    y     = as.vector(t(overlay_mat[plants_subset, , drop = FALSE]))
  )
  ov$plant <- factor(ov$plant, levels = plants_subset,
                     labels = plant_labels[plants_subset])
  ov$date  <- week_dates[ov$week]

  mm <- data.frame(
    plant = rep(plants_subset, each = T_wks),
    week  = rep(seq_len(T_wks), times = length(plants_subset)),
    miss  = as.vector(t(miss_mask[plants_subset, , drop = FALSE]))
  )
  mm <- mm[mm$miss, ]
  mm$plant <- factor(mm$plant, levels = plants_subset,
                     labels = plant_labels[plants_subset])
  mm$date  <- week_dates[mm$week]

  # Per-plant scaling of the overlay
  scale_df <- df |>
    group_by(plant) |>
    summarise(max_post = max(hi, na.rm = TRUE), .groups = "drop") |>
    left_join(
      ov |> group_by(plant) |>
        summarise(max_ovl = max(y, na.rm = TRUE), .groups = "drop"),
      by = "plant"
    ) |>
    mutate(scl = max_post / max_ovl)

  ov <- ov |>
    left_join(scale_df |> select(plant, scl), by = "plant") |>
    mutate(y_scl = y * scl)

  facet_levels <- levels(df$plant)
  scales_list <- lapply(facet_levels, function(pl_lbl) {
    scl_p <- scale_df$scl[scale_df$plant == pl_lbl]
    scale_y_continuous(sec.axis = sec_axis(~./scl_p, name = right_lab))
  })

  p <- ggplot(df, aes(x = date)) +
    geom_rect(data = mm,
              aes(xmin = date - 3.5, xmax = date + 3.5,
                  ymin = -Inf, ymax = Inf),
              fill = "grey85", alpha = 0.5, inherit.aes = FALSE) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = ribbon_color, alpha = 0.3) +
    geom_line(aes(y = md), color = ribbon_color, linewidth = 0.6) +
    geom_line(data = ov, aes(y = y_scl), color = "black",
              linewidth = 0.4, alpha = 0.8)

  if (add_hline_1) {
    p <- p + geom_hline(yintercept = 1, linetype = "dashed",
                        color = "grey30", linewidth = 0.7)
  }

  p +
    facet_wrap(~ plant, scales = "free_y", ncol = 5) +
    ggh4x::facetted_pos_scales(y = scales_list) +
    scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
    labs(x = NULL, y = y_lab) +
    theme_minimal(base_size = 10) +
    theme(strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

# ---- 4. Generate the 4 panel plots ----
cat("Generating RW2 panel plots...\n")

# Filter (ssm) variants
p_Rit_filter <- make_plot(Rt_ssm, y_filt, plant_names,
                          bquote(R[it]), "Filter (B gc/day)",
                          ribbon_color = COL_RT,
                          add_hline_1 = TRUE)
p_Iit_filter <- make_plot(It_ssm, y_filt, plant_names,
                          bquote(I[it]), "Filter (B gc/day)",
                          ribbon_color = COL_IT)

# Data (direct) variants
p_Rit_data <- make_plot(Rt_direct, y_raw, plant_names,
                        bquote(R[it]), "Viral Load (B gc/day)",
                        ribbon_color = COL_RT,
                        add_hline_1 = TRUE)
p_Iit_data <- make_plot(It_direct, y_raw, plant_names,
                        bquote(I[it]), "Viral Load (B gc/day)",
                        ribbon_color = COL_IT)

# ---- 5. Save with the RW2 naming convention ----
ggsave(file.path(OUT_DIR, "Rit_filter_rho_est_rw2.pdf"),
       p_Rit_filter, width = 15, height = 8)
ggsave(file.path(OUT_DIR, "Iit_filter_rho_est_rw2.pdf"),
       p_Iit_filter, width = 15, height = 8)
ggsave(file.path(OUT_DIR, "Rit_data_rho_est_rw2.pdf"),
       p_Rit_data,   width = 15, height = 8)
ggsave(file.path(OUT_DIR, "Iit_data_rho_est_rw2.pdf"),
       p_Iit_data,   width = 15, height = 8)
cat("  4 panel plots ->", OUT_DIR, "\n")

cat("\nDone.\n")
