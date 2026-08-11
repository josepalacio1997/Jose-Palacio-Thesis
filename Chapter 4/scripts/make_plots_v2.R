# =============================================================
# Modified panels + WWTP1 standalone plots
#   - fig7/8 direct_est (Data):  overlay = raw y_it
#   - fig7/8 ssm_est (Filter):   overlay = y_filt (Stage-1 smoothed)
#   - Gray bars for missing raw y in ALL plots
#   - Standalone WWTP1 plots for the 4 specifications
#
# Note: the .rds files are save_object lists, not cmdstanr fits.
# They contain: draws (draws_array), y_raw [P x T], y_filt [P x T].
# =============================================================

suppressPackageStartupMessages({
  library(posterior)
  library(ggplot2); library(dplyr); library(tidyr)
  library(ggh4x)      # per-facet secondary axes
  library(patchwork)  # combining plots
})

# Colors: Rt in orange, It in blue
COL_RT <- "darkorange"
COL_IT <- "steelblue"

# ---- PATHS ----
BASE <- "~/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW/Modelos AR con rho est. y rho=0"
FIT_SSM_EST     <- file.path(BASE, "fits/A_ssm_est.rds")
FIT_DIRECT_EST  <- file.path(BASE, "fits/A_direct_est.rds")
FIT_SSM_FIX     <- file.path(BASE, "fits/A_ssm_fix.rds")
FIT_DIRECT_FIX  <- file.path(BASE, "fits/A_direct_fix.rds")
OUT_DIR         <- file.path(BASE, "results_modelA_AR1")
POP_CACHE       <- file.path(BASE, "../outputs/N_plantas.rds")

# ---- 1. Load save_objects (contain draws + y_raw + y_filt) ----
obj_ssm        <- readRDS(FIT_SSM_EST)
obj_direct     <- readRDS(FIT_DIRECT_EST)
obj_ssm_fix    <- readRDS(FIT_SSM_FIX)
obj_direct_fix <- readRDS(FIT_DIRECT_FIX)
cat("Loaded fits.\n")

# y_raw / y_filt are the same in both variants (same underlying data);
# take them from obj_ssm because ssm variant needs y_filt.
y_raw  <- obj_ssm$y_raw    # [P x T]
y_filt <- obj_ssm$y_filt   # [P x T]

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

# ---- 1b. Build plant labels with populations "WWTP1 (578K)" style ----
plant_names <- paste0("WWTP", seq_len(m))
if (file.exists(POP_CACHE)) {
  N_plantas <- readRDS(POP_CACHE)
  plant_labels <- setNames(
    ifelse(is.na(N_plantas),
           plant_names,
           sprintf("%s (%dK)", plant_names, round(N_plantas / 1000))),
    plant_names
  )
  cat("Plant labels loaded with population sizes.\n")
} else {
  warning("N_plantas.rds not found -- falling back to bare WWTP labels")
  plant_labels <- setNames(plant_names, plant_names)
}
cat("Plant labels:\n"); print(plant_labels); cat("\n")

# ---- 2. Extract posterior summaries for Rt and It ----
# The draws are a draws_array with variables named "Rt[1,1]" etc.
extract_post <- function(save_obj, param) {
  d <- as_draws_matrix(save_obj$draws)
  vars <- grep(paste0("^", param, "\\["), colnames(d), value = TRUE)
  q <- apply(d[, vars, drop = FALSE], 2,
             quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  out <- data.frame(var = vars,
                    lo = q[1, ], md = q[2, ], hi = q[3, ])
  # Rt[p, t]  parse
  parsed <- do.call(rbind,
    lapply(strsplit(sub(".*\\[(.*)\\]", "\\1", vars), ","),
           function(x) as.integer(x)))
  out$plant_idx <- parsed[, 1]
  out$week      <- parsed[, 2]
  out$plant     <- plant_names[out$plant_idx]  # short name for indexing
  out
}

cat("Extracting Rt and It posteriors...\n")
Rt_ssm        <- extract_post(obj_ssm,        "Rt")
Rt_direct     <- extract_post(obj_direct,     "Rt")
Rt_ssm_fix    <- extract_post(obj_ssm_fix,    "Rt")
Rt_direct_fix <- extract_post(obj_direct_fix, "Rt")
It_ssm        <- extract_post(obj_ssm,        "It")
It_direct     <- extract_post(obj_direct,     "It")
It_ssm_fix    <- extract_post(obj_ssm_fix,    "It")
It_direct_fix <- extract_post(obj_direct_fix, "It")
cat("  done.\n")

# ---- 3. Plot function ----
# plants_subset uses short names ("WWTP1", ...). Facet labels are the
# enriched labels ("WWTP1 (578K)", ...) via plant_labels[].
# Per-plant scaling of the overlay so each panel is readable.
# If y_lim is provided (e.g., c(0, 3)), the left y-axis is fixed to that
# range for every panel, and the overlay is scaled to fit within it.
make_plot <- function(post_df, overlay_mat, plants_subset,
                      y_lab, right_lab, ribbon_color = "steelblue",
                      title_extra = "", y_lim = NULL,
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

  # ---- Per-plant scaling ----
  # If y_lim is set, use y_lim[2] as the reference; else use the panel's own hi max.
  scale_df <- df |>
    group_by(plant) |>
    summarise(max_post = max(hi, na.rm = TRUE), .groups = "drop") |>
    left_join(
      ov |> group_by(plant) |>
        summarise(max_ovl = max(y, na.rm = TRUE), .groups = "drop"),
      by = "plant"
    ) |>
    mutate(scl = if (is.null(y_lim)) max_post / max_ovl
                 else y_lim[2] / max_ovl)

  ov <- ov |>
    left_join(scale_df |> select(plant, scl), by = "plant") |>
    mutate(y_scl = y * scl)

  # Per-facet secondary axes (one scale object per plant)
  facet_levels <- levels(df$plant)
  scales_list <- lapply(facet_levels, function(pl_lbl) {
    scl_p <- scale_df$scl[scale_df$plant == pl_lbl]
    if (is.null(y_lim)) {
      scale_y_continuous(sec.axis = sec_axis(~./scl_p, name = right_lab))
    } else {
      scale_y_continuous(limits = y_lim,
                         sec.axis = sec_axis(~./scl_p, name = right_lab))
    }
  })

  # If y_lim is provided, use "fixed" scales; otherwise free_y (per-panel)
  facet_scales <- if (is.null(y_lim)) "free_y" else "fixed"

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
    facet_wrap(~ plant, scales = facet_scales, ncol = 5) +
    ggh4x::facetted_pos_scales(y = scales_list) +
    scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
    labs(x = NULL, y = y_lab, title = title_extra) +
    theme_minimal(base_size = 10) +
    theme(strip.text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

# ---- 4. Panel plots (15 WWTPs each) ----
# 8 plots total: 4 variants (Filter/Data x est/fix) x 2 params (Rt/It).
# Rt plots include a dashed horizontal line at R = 1.
# Filter variants overlay y_filt; Data variants overlay y_raw.
cat("Generating panel plots...\n")

# ---- Filter variants (overlay = y_filt) ----
p_Rit_filter_rho_est <- make_plot(Rt_ssm, y_filt, plant_names,
                                  bquote(R[it]), "Filter (B gc/day)",
                                  ribbon_color = COL_RT,
                                  add_hline_1 = TRUE)
p_Iit_filter_rho_est <- make_plot(It_ssm, y_filt, plant_names,
                                  bquote(I[it]), "Filter (B gc/day)",
                                  ribbon_color = COL_IT)
p_Rit_filter_rho_0   <- make_plot(Rt_ssm_fix, y_filt, plant_names,
                                  bquote(R[it]), "Filter (B gc/day)",
                                  ribbon_color = COL_RT,
                                  add_hline_1 = TRUE)
p_Iit_filter_rho_0   <- make_plot(It_ssm_fix, y_filt, plant_names,
                                  bquote(I[it]), "Filter (B gc/day)",
                                  ribbon_color = COL_IT)

# ---- Data variants (overlay = y_raw) ----
p_Rit_data_rho_est <- make_plot(Rt_direct, y_raw, plant_names,
                                bquote(R[it]), "Viral Load (B gc/day)",
                                ribbon_color = COL_RT,
                                add_hline_1 = TRUE)
p_Iit_data_rho_est <- make_plot(It_direct, y_raw, plant_names,
                                bquote(I[it]), "Viral Load (B gc/day)",
                                ribbon_color = COL_IT)
p_Rit_data_rho_0   <- make_plot(Rt_direct_fix, y_raw, plant_names,
                                bquote(R[it]), "Viral Load (B gc/day)",
                                ribbon_color = COL_RT,
                                add_hline_1 = TRUE)
p_Iit_data_rho_0   <- make_plot(It_direct_fix, y_raw, plant_names,
                                bquote(I[it]), "Viral Load (B gc/day)",
                                ribbon_color = COL_IT)

# ---- Save each as its own PDF with the short naming convention ----
ggsave(file.path(OUT_DIR, "Rit_filter_rho_est.pdf"),
       p_Rit_filter_rho_est, width = 15, height = 8)
ggsave(file.path(OUT_DIR, "Iit_filter_rho_est.pdf"),
       p_Iit_filter_rho_est, width = 15, height = 8)
ggsave(file.path(OUT_DIR, "Rit_filter_rho_0.pdf"),
       p_Rit_filter_rho_0,   width = 15, height = 8)
ggsave(file.path(OUT_DIR, "Iit_filter_rho_0.pdf"),
       p_Iit_filter_rho_0,   width = 15, height = 8)
ggsave(file.path(OUT_DIR, "Rit_data_rho_est.pdf"),
       p_Rit_data_rho_est,   width = 15, height = 8)
ggsave(file.path(OUT_DIR, "Iit_data_rho_est.pdf"),
       p_Iit_data_rho_est,   width = 15, height = 8)
ggsave(file.path(OUT_DIR, "Rit_data_rho_0.pdf"),
       p_Rit_data_rho_0,     width = 15, height = 8)
ggsave(file.path(OUT_DIR, "Iit_data_rho_0.pdf"),
       p_Iit_data_rho_0,     width = 15, height = 8)
cat("  8 panel plots ->", OUT_DIR, "\n")

# ---- 5. WWTP1 combined 2x1 plots ----
# Top row = Viral Load overlay (direct_est / raw data)
# Bottom row = Filter overlay (ssm_est / MARSS)
cat("Generating WWTP1 combined plots...\n")

# Rt combined -- left y-axis fixed to [0, 3] on both panels
p_rt_top <- make_plot(Rt_direct, y_raw,  "WWTP1", bquote(R[it]),
                      "Viral Load (B gc/day)",
                      ribbon_color = COL_RT, title_extra = NULL,
                      y_lim = c(0, 4))
p_rt_bot <- make_plot(Rt_ssm,    y_filt, "WWTP1", bquote(R[it]),
                      "Filter (B gc/day)",
                      ribbon_color = COL_RT, title_extra = NULL,
                      y_lim = c(0, 4))
combined_rt <- p_rt_top / p_rt_bot
ggsave(file.path(OUT_DIR, "WWTP1_Rt_combined.pdf"),
       combined_rt, width = 6, height = 8)
cat("  -> WWTP1_Rt_combined.pdf\n")

# It combined
p_it_top <- make_plot(It_direct, y_raw,  "WWTP1", bquote(I[it]),
                      "Viral Load (B gc/day)",
                      ribbon_color = COL_IT, title_extra = NULL)
p_it_bot <- make_plot(It_ssm,    y_filt, "WWTP1", bquote(I[it]),
                      "Filter (B gc/day)",
                      ribbon_color = COL_IT, title_extra = NULL)
combined_it <- p_it_top / p_it_bot
ggsave(file.path(OUT_DIR, "WWTP1_It_combined.pdf"),
       combined_it, width = 6, height = 8)
cat("  -> WWTP1_It_combined.pdf\n")

# ---- 6. WWTP1 standalone for Filter with rho = 0 ----
cat("Generating WWTP1 standalone plots for Filter with rho = 0...\n")
p_rt_ssm_fix <- make_plot(Rt_ssm_fix, y_filt, "WWTP1", bquote(R[it]),
                          "Filter (B gc/day)",
                          ribbon_color = COL_RT, title_extra = NULL,
                          y_lim = c(0, 4))
ggsave(file.path(OUT_DIR, "WWTP1_Rt_ssm_fix.pdf"),
       p_rt_ssm_fix, width = 6, height = 4)
cat("  -> WWTP1_Rt_ssm_fix.pdf\n")

p_it_ssm_fix <- make_plot(It_ssm_fix, y_filt, "WWTP1", bquote(I[it]),
                          "Filter (B gc/day)",
                          ribbon_color = COL_IT, title_extra = NULL)
ggsave(file.path(OUT_DIR, "WWTP1_It_ssm_fix.pdf"),
       p_it_ssm_fix, width = 6, height = 4)
cat("  -> WWTP1_It_ssm_fix.pdf\n")

cat("\nDone. All outputs in: ", OUT_DIR, "\n")
