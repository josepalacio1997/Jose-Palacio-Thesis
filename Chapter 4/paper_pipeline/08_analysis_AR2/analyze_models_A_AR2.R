# =============================================================================
# analyze_models_A_AR2.R
# -----------------------------------------------------------------------------
# second-order-only analysis. Mirror of analyze_models_A_AR1.R but for the
# second-order variants of the latent random walk on lambda. Writes outputs to a
# DEDICATED folder (results_modelA_AR2/) so first-order and second-order results never mix.
#
# Expected fits (4 second-order variants):
#   A_ssm_est_ar2   : SSM Stage-1 filter, phi estimado (hat(rho))
#   A_direct_est_ar2 : Direct, hat(rho) (hat(rho))
#   A_ssm_fix_ar2   : SSM, phi = 0       (may be MISSING if not yet produced)
#   A_direct_fix_ar2 : Direct, rho = 0     (may be MISSING if not yet produced)
#
# Note: comparison against EpiSewer (single-plant baseline) is NOT replicated
# here because the headline paper comparison is against the first-order Direct
# variant. second-order results live in the appendix as a sensitivity check.
#
# Output figures and tables (in results_modelA_AR2/):
#   fig6_phi_posterior_AR2.pdf       - density of phi posterior, second-order
#   fig7_Rt_by_plant_AR2.pdf         - Rt trajectories per plant, second-order
#   fig8_It_by_plant_AR2.pdf         - It trajectories per plant, second-order
#   fig11_rhat_by_block_AR2.pdf      - Rhat diagnostics
#   fig12_ess_by_block_AR2.pdf       - ESS diagnostics
#   table1_beta_posterior_AR2.csv    - per-plant beta posterior
#   tableS_phi_posterior_AR2.csv     - phi posterior summary
#   convergence_diagnostic_AR2.csv   - divergences, treedepth, ebfmi
# =============================================================================

suppressPackageStartupMessages({
  library(posterior); library(loo)
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(ggplot2); library(scales)
  library(grid); library(gridExtra)
  # ggh4x para sec_axis por panel con free_y
  if (!requireNamespace("ggh4x", quietly = TRUE)) {
    install.packages("ggh4x", repos = "https://cloud.r-project.org")
  }
  library(ggh4x)
  # scoringRules para CRPS y otras proper scoring rules
  if (!requireNamespace("scoringRules", quietly = TRUE)) {
    install.packages("scoringRules", repos = "https://cloud.r-project.org")
  }
  library(scoringRules)
})

PROJ_DIR <- "/Users/josepalacio/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW"
OUT_DIR  <- file.path(PROJ_DIR, "Modelos AR con rho est. y rho=0",
                      "results_modelA_AR2")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# House style (Rice + HHD colors)
source(file.path(PROJ_DIR, "Modelos AR con rho est. y rho=0", "R",
                 "house_style.R"))

# ============================================================================
# 0 . Plant populations + labels (igual que analyze_stan_results.qmd del paper)
# ============================================================================
P <- 15
plant_names <- paste0("WWTP", 1:P)

pop_cache <- file.path(PROJ_DIR, "outputs", "N_plantas.rds")
if (file.exists(pop_cache)) {
  N_plantas <- readRDS(pop_cache)
  cat("Poblaciones cargadas de cache N_plantas.rds\n")
} else {
  stop("N_plantas.rds no existe. Corre primero analyze_stan_results.qmd para generar.")
}

# Etiquetas "WWTP1 (578K)" estilo paper
plant_labels <- setNames(
  ifelse(is.na(N_plantas),
         plant_names,
         sprintf("%s (%dK)", plant_names, round(N_plantas / 1000))),
  plant_names
)
cat("Plant labels:\n"); print(plant_labels); cat("\n")

# ============================================================================
# 1 . Load fits (6 fits: first-order x SSM/Direct x rho-est/rho-fix +
#                                 second-order x SSM/Direct x rho-est solamente)
#     Symlinks en Modelos AR.../fits/
#
#     SSM   = Stage-1 state-space filter (Kalman) supplies sigma_obs externally
#     Direct = sigma_obs estimated jointly with the latent process inside Stan
#
#     Nota: second-order rho=0 fue descartado (no se corrio en Mac por geometria lenta).
# ============================================================================
FITS_DIR <- file.path(PROJ_DIR, "Modelos AR con rho est. y rho=0", "fits")
FIT_PATHS <- list(
  # ---- second-order only ----
  A_ssm_est_ar2   = file.path(FITS_DIR, "A_ssm_est_ar2.rds"),
  A_direct_est_ar2 = file.path(FITS_DIR, "A_direct_est_ar2.rds"),
  A_ssm_fix_ar2   = file.path(FITS_DIR, "A_ssm_fix_ar2.rds"),
  A_direct_fix_ar2 = file.path(FITS_DIR, "A_direct_fix_ar2.rds")
)

cat("============================================================\n")
cat("  Loading model fits
")
cat("============================================================\n\n")
fits <- list()
for (lab in names(FIT_PATHS)) {
  p <- FIT_PATHS[[lab]]
  if (file.exists(p)) {
    fits[[lab]] <- readRDS(p)
    sz <- file.size(p) / 1e6
    cat(sprintf("  OK  %-15s  [%.0f MB]\n", lab, sz))
  } else {
    cat(sprintf("  --  %-15s  MISSING: %s\n", lab, p))
  }
}
cat(sprintf("\n  Loaded %d / %d fits\n\n", length(fits), length(FIT_PATHS)))
stopifnot(length(fits) >= 1)

# Colores second-order: orange for SSM, blue for Direct; tono claro = rho-fix, oscuro = rho-est
COL_SSM_EST_AR2   <- HHD_ORANGE        # SSM rho-est
COL_DIRECT_EST_AR2 <- "#2C7BB6"          # Direct rho-est: azul medio
COL_SSM_FIX_AR2   <- HHD_ORANGE_LIGHT   # SSM rho=0
COL_DIRECT_FIX_AR2 <- "#9ECAE1"          # Direct rho=0: azul claro

COLORS_4 <- c(
  A_ssm_est_ar2   = COL_SSM_EST_AR2,
  A_direct_est_ar2 = COL_DIRECT_EST_AR2,
  A_ssm_fix_ar2   = COL_SSM_FIX_AR2,
  A_direct_fix_ar2 = COL_DIRECT_FIX_AR2
)

# Titulos LaTeX via bquote(); rho-estimated -> hat(rho), rho-fix -> rho == 0
title_label <- list(
  A_ssm_est_ar2   = bquote("Filter; " * hat(rho)),
  A_direct_est_ar2 = bquote("Data; " * hat(rho)),
  A_ssm_fix_ar2   = bquote("Filter; " * rho == 0),
  A_direct_fix_ar2 = bquote("Data; " * rho == 0)
)

# Version texto plano (CSV, consola)
title_label_plain <- c(
  A_ssm_est_ar2   = "Filter; rho-est",
  A_direct_est_ar2 = "Data; rho-est",
  A_ssm_fix_ar2   = "Filter; rho=0",
  A_direct_fix_ar2 = "Data; rho=0"
)
# Plotmath strings for axis ticks / legend keys (use with parse(text=...))
title_label_parse <- c(
  A_ssm_est_ar2    = "'Filter; '*hat(rho)",
  A_direct_est_ar2 = "'Data; '*hat(rho)",
  A_ssm_fix_ar2    = "'Filter; '*rho==0",
  A_direct_fix_ar2 = "'Data; '*rho==0"
)
# Helper: apply parse(text=...) lookup to a character vector of fit IDs
fit_label_parser <- function(x) {
  parse(text = title_label_parse[as.character(x)])
}

# Plotmath labeller for the parameter block facet (beta, It, Rt) so the
# strip text in fig11/fig12/fig13 renders as math (beta[i], I[it], R[it])
# instead of plain "beta", "It", "Rt".
block_label_parse <- c(
  "beta" = "beta[i]",
  "It"   = "I[it]",
  "Rt"   = "R[it]",
  "rho"  = "rho"
)
block_labeller <- as_labeller(block_label_parse, label_parsed)

# ============================================================================
# 2 . Helpers para extraer drawsy quantiles posteriores
# ============================================================================
get_draws_matrix <- function(fit, vars) {
  d <- fit$draws
  if (!is.null(vars)) d <- posterior::subset_draws(d, variable = vars)
  posterior::as_draws_matrix(d)
}

# Para Rt[p, t] o It[p, t]: extraer draws y armar long df con quantiles 50%+95%
extract_PT_posterior <- function(fit, varname, fit_label) {
  dm <- get_draws_matrix(fit, varname)
  df <- data.frame(
    variable = colnames(dm),
    mean     = colMeans(dm),
    q025     = apply(dm, 2, quantile, probs = 0.025, na.rm = TRUE),
    q25      = apply(dm, 2, quantile, probs = 0.25,  na.rm = TRUE),
    q50      = apply(dm, 2, quantile, probs = 0.50,  na.rm = TRUE),
    q75      = apply(dm, 2, quantile, probs = 0.75,  na.rm = TRUE),
    q975     = apply(dm, 2, quantile, probs = 0.975, na.rm = TRUE)
  )
  pat <- sprintf("%s\\[(\\d+),(\\d+)\\]", varname)
  idx <- str_match(df$variable, pat)
  df$p <- as.integer(idx[, 2])
  df$t <- as.integer(idx[, 3])
  df$plant <- factor(plant_names[df$p],
                     levels = plant_names,
                     labels = plant_labels[plant_names])
  df$fit <- fit_label
  df
}

# Per-plant beta posterior
extract_beta_posterior <- function(fit, fit_label) {
  bm <- get_draws_matrix(fit, "beta")
  P_loc <- ncol(bm)
  data.frame(
    fit       = fit_label,
    plant     = plant_names[seq_len(P_loc)],
    plant_lab = plant_labels[seq_len(P_loc)],
    beta_mean = colMeans(bm),
    beta_q025 = apply(bm, 2, quantile, probs = 0.025, na.rm = TRUE),
    beta_med  = apply(bm, 2, quantile, probs = 0.50,  na.rm = TRUE),
    beta_q975 = apply(bm, 2, quantile, probs = 0.975, na.rm = TRUE),
    row.names = NULL
  )
}

# ============================================================================
# 3 . Convergence diagnostics
# ============================================================================
cat("== 3 . Convergencia ==\n")
diag_rows <- list()
for (lab in names(fits)) {
  d <- fits[[lab]]$diagnostic_summary
  meta <- fits[[lab]]$meta
  diag_rows[[lab]] <- data.frame(
    fit            = lab,
    iter_sampling  = if (!is.null(meta$iter_sampling)) meta$iter_sampling else NA,
    num_divergent  = sum(d$num_divergent),
    num_treedepth  = sum(d$num_max_treedepth),
    ebfmi_min      = if (!is.null(d$ebfmi)) round(min(d$ebfmi), 3) else NA
  )
}
diag_tbl <- do.call(rbind, diag_rows); rownames(diag_tbl) <- NULL
print(diag_tbl, row.names = FALSE)
write.csv(diag_tbl, file.path(OUT_DIR, "convergence_diagnostic_AR2.csv"),
          row.names = FALSE)
cat("\n")

# ============================================================================
# 4 . Phi posterior (solo rho-est fits) + Figure 6 (densidad overlay)
# ============================================================================
cat("== 4 . Phi posterior + Figure 6 ==\n")

phi_summ <- list(); phi_dens <- list()
for (lab in names(fits)) {
  r <- fits[[lab]]$rho_posterior_summary
  if (!is.null(r)) {
    phi_summ[[lab]] <- data.frame(
      fit  = lab,
      mean = r$mean, sd = r$sd,
      q5   = r$`5%`, q50 = r$`50%`, q95 = r$`95%`,
      rhat = r$rhat, ess_bulk = r$ess_bulk
    )
    rho_draws <- as.numeric(get_draws_matrix(fits[[lab]], "rho"))
    phi_dens[[lab]] <- data.frame(fit = lab, phi = rho_draws)
  }
}
phi_tbl <- do.call(rbind, phi_summ); rownames(phi_tbl) <- NULL
print(phi_tbl, row.names = FALSE)
write.csv(phi_tbl, file.path(OUT_DIR, "tableS_phi_posterior_AR2.csv"),
          row.names = FALSE)

if (length(phi_dens) > 0) {
  phi_df <- bind_rows(phi_dens)
  est_fits <- intersect(c("A_ssm_est_ar2", "A_direct_est_ar2"),
                         unique(phi_df$fit))
  # Keep fit IDs; the legend will parse them via title_label_parse.
  phi_df$fit <- factor(phi_df$fit, levels = est_fits)
  phi_means <- phi_df |> group_by(fit) |>
    summarise(mn = mean(phi), .groups = "drop")

  p_fig6 <- ggplot(phi_df, aes(x = phi, fill = fit, color = fit)) +
    geom_density(alpha = 0.35, linewidth = 0.7) +
    geom_vline(data = phi_means, aes(xintercept = mn, color = fit),
               linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
    scale_fill_manual(values = COLORS_4[est_fits], name = NULL,
                      labels = fit_label_parser) +
    scale_color_manual(values = COLORS_4[est_fits], name = NULL,
                       labels = fit_label_parser) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(
      title    = NULL,
      subtitle = NULL,
      x = expression(rho), y = "Posterior density"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title    = element_text(face = "bold", color = RICE_BLUE_DARK),
          plot.subtitle = element_text(color = NEUTRAL_GRAY))

  ggsave(file.path(OUT_DIR, "fig6_phi_posterior_AR2.pdf"),
         p_fig6, width = 7, height = 4.5)
  cat("  -> fig6_phi_posterior.pdf\n")
}
cat("\n")

# ============================================================================
# 5 . Viral load + missing data extraction (real y_raw, independiente del modelo)
# ============================================================================
cat("== 5 . Extracting viral loads (y_raw, B gc/day) y missing data mask ==\n")
# y_raw esta guardado en cualquier fit (matriz P x T). Tomamos del primer fit.
y_raw <- fits[[1]]$y_raw
stopifnot(is.matrix(y_raw), nrow(y_raw) == P)
T_focus <- ncol(y_raw)
cat(sprintf("  y_raw: %d plants x %d weeks; %d NAs (%.1f%%)\n",
            P, T_focus, sum(is.na(y_raw)),
            100 * mean(is.na(y_raw))))

# Long format
viral_load_long <- as.data.frame(y_raw) |>
  setNames(seq_len(T_focus)) |>
  mutate(p = seq_len(P)) |>
  pivot_longer(-p, names_to = "t", values_to = "viral_load") |>
  mutate(
    t          = as.integer(t),
    plant      = factor(plant_names[p],
                        levels = plant_names,
                        labels = plant_labels[plant_names]),
    is_missing = is.na(viral_load)
  )

# Missing data df: una fila por (plant, t) donde y_raw es NA
missing_long <- viral_load_long |> filter(is_missing)
cat(sprintf("  Missing rows: %d (de %d totales)\n",
            nrow(missing_long), nrow(viral_load_long)))

# ============================================================================
# 6 . Figure 7 : Rt by plant (estilo paper, con viral load y missing bars)
# ============================================================================
cat("== 6 . Figure 7 : Rt by plant ==\n")
rt_long <- imap_dfr(fits, function(fit, lab) {
  cat(sprintf("  Extracting Rt para %s ...\n", lab))
  extract_PT_posterior(fit, "Rt", lab)
})
write.csv(rt_long |> select(fit, plant, p, t, mean, q025, q25, q50, q75, q975),
          file.path(OUT_DIR, "Rt_posterior_AR2.csv"), row.names = FALSE)

# ----------------------------------------------------------------------------
# make_paper_panel: panel multi-planta estilo paper
#   - posterior median + bandas 50%/95%
#   - viral load real (y_raw) escalada por planta, linea negra solida
#   - barras grises verticales en weeks con y_raw = NA
#   - sec_axis por panel via ggh4x para mostrar valores reales en B gc/day
# ----------------------------------------------------------------------------
make_paper_panel <- function(df, color_fill, title, y_lab,
                              add_hline_1 = FALSE) {

  # 1. Per-plant scaling: viral_load * scale = mismo rango que q975
  vl_summary <- viral_load_long |>
    group_by(plant) |>
    summarise(
      vl_max = suppressWarnings(max(viral_load, na.rm = TRUE)),
      .groups = "drop"
    )
  yt_summary <- df |>
    group_by(plant) |>
    summarise(
      yt_max = suppressWarnings(max(q975, na.rm = TRUE)),
      yt_min = suppressWarnings(min(q025, na.rm = TRUE)),
      .groups = "drop"
    )
  scale_per_plant <- yt_summary |>
    left_join(vl_summary, by = "plant") |>
    mutate(
      scale = ifelse(is.finite(vl_max) & vl_max > 0,
                      yt_max / vl_max, 1),
      scale = ifelse(is.finite(scale) & scale > 0, scale, 1)
    )

  # 2. Viral load scaled to fit each panel's range
  vl_for_plot <- viral_load_long |>
    left_join(scale_per_plant |> select(plant, scale, vl_max),
              by = "plant") |>
    mutate(vl_scaled = viral_load * scale)

  # 3. Missing data: per-plant ymin/ymax para que cubra todo el panel
  missing_for_plot <- missing_long |>
    left_join(yt_summary, by = "plant") |>
    mutate(
      ymin = ifelse(is.finite(yt_min), yt_min, -Inf),
      ymax = ifelse(is.finite(yt_max), yt_max * 1.1, Inf)
    )

  p <- ggplot(df, aes(x = t)) +
    # Barras grises verticales para missing data (atras de todo)
    geom_rect(data = missing_for_plot,
              aes(xmin = t - 0.5, xmax = t + 0.5,
                  ymin = -Inf,    ymax = Inf),
              fill = "grey85", alpha = 0.45,
              inherit.aes = FALSE) +
    # Bandas 95% (claro) y 50% (oscuro) del posterior
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = color_fill,
                alpha = 0.18, color = NA) +
    geom_ribbon(aes(ymin = q25,  ymax = q75),  fill = color_fill,
                alpha = 0.40, color = NA) +
    geom_line(aes(y = q50), color = color_fill, linewidth = 0.55) +
    # Viral load real, escalado per-plant, linea negra solida
    geom_line(data = vl_for_plot,
              aes(x = t, y = vl_scaled),
              color = "black", linewidth = 0.40,
              inherit.aes = FALSE) +
    facet_wrap(~ plant, ncol = 3, scales = "free_y") +
    labs(x = "Week", y = y_lab) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none",
          strip.text     = element_text(face = "bold",
                                          color = RICE_BLUE_DARK),
          plot.title     = element_text(face = "bold",
                                          color = RICE_BLUE_DARK),
          panel.spacing  = unit(0.5, "lines"))

  if (add_hline_1) {
    p <- p + geom_hline(yintercept = 1, linetype = "dashed",
                          color = "grey30", linewidth = 0.7)
  }

  # 4. sec_axis por panel via ggh4x::facetted_pos_scales().
  # Por cada planta, transform = ~ . / scale_factor[planta] -> recupera
  # los valores reales de viral load en B gc/day.
  plant_levels <- levels(df$plant)
  scale_lookup <- setNames(scale_per_plant$scale, scale_per_plant$plant)

  facet_scales <- lapply(plant_levels, function(pl) {
    sc <- scale_lookup[[pl]]
    if (!is.finite(sc) || sc <= 0) sc <- 1
    scale_y_continuous(
      sec.axis = sec_axis(
        transform = ~ . / sc,
        name      = "Viral load (B gc/day)",
        breaks    = scales::pretty_breaks(n = 4)
      )
    )
  })

  p + facetted_pos_scales(y = facet_scales)
}

# ----------------------------------------------------------------------------
# make_multipanel_pdf: una pagina por fit dentro del PDF
# ----------------------------------------------------------------------------
# Color convention: single tone per parameter across all variants.
# Rt -> orange (same shade for est and fix); It -> blue (same shade).
pick_color <- function(varname, lab) {
  if (varname == "Rt") HHD_ORANGE else "#2C7BB6"
}

make_multipanel_pdf <- function(long_df, varname, file_out, y_lab) {
  pdf(file_out, width = 16, height = 12)
  for (lab in names(fits)) {
    sub <- long_df |> filter(fit == lab)
    fit_tag <- title_label[[lab]]
    base_plot <- make_paper_panel(
      sub, pick_color(varname, lab),
      title       = fit_tag,
      y_lab       = y_lab,
      add_hline_1 = (varname == "Rt")
    )
    print(base_plot)
  }
  dev.off()
}

# Helper: emit ONE single-page PDF per variant, for inclusion in the appendix
# as separate figures. Filenames follow the pattern <prefix>_<variant>.pdf
make_per_variant_pdfs <- function(long_df, varname, prefix, y_lab) {
  paths <- character(0)
  for (lab in names(fits)) {
    sub <- long_df |> filter(fit == lab)
    fit_tag <- title_label[[lab]]
    base_plot <- make_paper_panel(
      sub, pick_color(varname, lab),
      title       = fit_tag,
      y_lab       = y_lab,
      add_hline_1 = (varname == "Rt")
    )
    out_path <- file.path(OUT_DIR, sprintf("%s_%s.pdf", prefix, lab))
    ggsave(out_path, base_plot, width = 16, height = 12, device = "pdf")
    paths <- c(paths, basename(out_path))
  }
  paths
}

make_multipanel_pdf(rt_long, "Rt",
                     file.path(OUT_DIR, "fig7_Rt_by_plant_AR2.pdf"),
                     y_lab = expression(R[t]))
cat("  -> fig7_Rt_by_plant.pdf (1 pagina por fit)\n")
rt_per_variant <- make_per_variant_pdfs(rt_long, "Rt",
                                          "fig7_Rt_by_plant_AR2",
                                          y_lab = expression(R[t]))
cat(sprintf("  -> per-variant: %s\n\n", paste(rt_per_variant, collapse = ", ")))

# ============================================================================
# 7 . Figure 8 : It by plant (1 pagina por fit, con viral load + missing bars)
# ============================================================================
cat("== 7 . Figure 8 : It by plant ==\n")
it_long <- imap_dfr(fits, function(fit, lab) {
  cat(sprintf("  Extracting It para %s ...\n", lab))
  extract_PT_posterior(fit, "It", lab)
})
write.csv(it_long |> select(fit, plant, p, t, mean, q025, q25, q50, q75, q975),
          file.path(OUT_DIR, "It_posterior_AR2.csv"), row.names = FALSE)

make_multipanel_pdf(it_long, "It",
                     file.path(OUT_DIR, "fig8_It_by_plant_AR2.pdf"),
                     y_lab = expression(I[t]))
cat("  -> fig8_It_by_plant.pdf (1 pagina por fit)\n")
it_per_variant <- make_per_variant_pdfs(it_long, "It",
                                          "fig8_It_by_plant_AR2",
                                          y_lab = expression(I[t]))
cat(sprintf("  -> per-variant: %s\n\n", paste(it_per_variant, collapse = ", ")))

# ============================================================================
# 8 . Tabla 1 : beta posterior mean + 95% CI por planta y paradigma
# ============================================================================
cat("== 8 . Tabla 1 : beta posterior ==\n")
beta_all <- imap_dfr(fits, function(fit, lab) extract_beta_posterior(fit, lab))

# Wide para tabla estilo paper: incluye los 4 fits rho-est (AR1 y AR2 x SSM y Joint)
est_fits_avail <- intersect(c("A_ssm_est", "A_direct_est",
                               "A_ssm_est_ar2", "A_direct_est_ar2"),
                             names(fits))
beta_wide <- beta_all |>
  filter(fit %in% est_fits_avail) |>
  mutate(
    param  = sprintf("%.2f", beta_mean),
    ci     = sprintf("[%.2f, %.2f]", beta_q025, beta_q975)
  ) |>
  select(plant_lab, fit, beta_mean, beta_q025, beta_q975) |>
  pivot_wider(
    names_from  = fit,
    values_from = c(beta_mean, beta_q025, beta_q975),
    names_glue  = "{fit}_{.value}"
  )

# Tambien una version completa con los 4 fits
write.csv(beta_all, file.path(OUT_DIR, "beta_posterior_long_AR2.csv"),
          row.names = FALSE)
write.csv(beta_wide, file.path(OUT_DIR, "table1_beta_posterior_AR2.csv"),
          row.names = FALSE)

print(beta_wide, row.names = FALSE, digits = 3)
cat("\n")
cat("  -> table1_beta_posterior.csv\n")
cat("  -> beta_posterior_long.csv\n\n")

# ============================================================================
# 8 . PSIS-LOO compare
# ============================================================================
cat("== 8 . PSIS-LOO compare ==\n")
loo_objs <- tryCatch(
  lapply(fits, function(fit) loo(get_draws_matrix(fit, "log_lik"))),
  error = function(e) {
    cat("  WARN: log_lik missing en alguno de los fits — skip LOO\n")
    NULL
  })
if (!is.null(loo_objs)) {
  loo_cmp  <- loo_compare(loo_objs)
  print(loo_cmp)
  loo_df <- as.data.frame(loo_cmp)
  loo_df$model <- rownames(loo_cmp)
  loo_df <- loo_df[, c("model", setdiff(names(loo_df), "model"))]
  write.csv(loo_df, file.path(OUT_DIR, "loo_compare_AR2.csv"), row.names = FALSE)
}
cat("\n")

# ============================================================================
# 9 . Correlacion entre todos los pares de fits (Rt y It scatter matrices)
#     Cada par responde una pregunta:
#       est vs fix (within paradigm) - ¿phi cambia R o I dentro del paradigma?
#       first-order vs second-order (cross AR) - ¿el orden del AR cambia las inferencias?
# ============================================================================
cat("== 9 . Agreement / correlation entre todos los pares ==\n")

# Posterior means por (fit, plant, t) — uso q50 para robustez
rt_mean <- rt_long |>
  select(fit, p, t, val = q50) |>
  pivot_wider(names_from = fit, values_from = val)
it_mean <- it_long |>
  select(fit, p, t, val = q50) |>
  pivot_wider(names_from = fit, values_from = val)

# Helper: pair statistics
pair_stats <- function(x, y, log_scale = FALSE) {
  if (log_scale) {
    keep <- is.finite(x) & is.finite(y) & x > 0 & y > 0
    lx <- log10(x[keep]); ly <- log10(y[keep])
    list(cor = cor(lx, ly), rmse = sqrt(mean((lx - ly)^2)),
         ratio = 10^median(ly - lx))
  } else {
    keep <- is.finite(x) & is.finite(y)
    xx <- x[keep]; yy <- y[keep]
    list(cor = cor(xx, yy), rmse = sqrt(mean((xx - yy)^2)),
         ratio = median(yy / pmax(xx, 1e-10)))
  }
}

combos <- combn(names(fits), 2, simplify = FALSE)

# Tabla de agreement
agree_rows <- list()
for (cb in combos) {
  a <- cb[1]; b <- cb[2]
  R_a <- pair_stats(rt_mean[[a]], rt_mean[[b]])
  I_a <- pair_stats(it_mean[[a]], it_mean[[b]], log_scale = TRUE)
  agree_rows[[paste(a, b, sep = "_vs_")]] <- data.frame(
    pair        = paste(a, "vs", b),
    cor_R       = round(R_a$cor, 4),
    rmse_R      = round(R_a$rmse, 4),
    ratio_R     = round(R_a$ratio, 4),
    cor_I_log10 = round(I_a$cor, 4),
    rmse_I_log  = round(I_a$rmse, 4),
    ratio_I     = round(I_a$ratio, 4)
  )
}
agree_tbl <- do.call(rbind, agree_rows); rownames(agree_tbl) <- NULL
print(agree_tbl, row.names = FALSE); cat("\n")
write.csv(agree_tbl, file.path(OUT_DIR, "agreement_all_pairs_AR2.csv"),
          row.names = FALSE)
cat("  -> agreement_all_pairs.csv\n")

# --- Figura: scatter Rt matrix (6 paneles, uno por par) ---
make_scatter_panel <- function(df, a, b, log_log = FALSE,
                                title_extra = "",
                                xy_lim = NULL) {
  d <- data.frame(x = df[[a]], y = df[[b]]) |>
    filter(is.finite(x) & is.finite(y))
  if (log_log) d <- filter(d, x > 0, y > 0)
  r <- if (log_log) cor(log10(d$x), log10(d$y)) else cor(d$x, d$y)

  # Map fit IDs to formatted axis labels (bquote -> rho-hat, rho==0 etc.)
  x_lab <- title_label[[a]]; if (is.null(x_lab)) x_lab <- a
  y_lab <- title_label[[b]]; if (is.null(y_lab)) y_lab <- b
  if (log_log) {
    xx <- log10(d$x); yy <- log10(d$y)
  } else {
    xx <- d$x; yy <- d$y
  }
  r   <- cor(xx, yy)
  ann <- sprintf("r = %.3f", r)
  # Top-left annotation: data-scale coords for log_log, raw for linear
  if (!is.null(xy_lim)) {
    if (log_log) {
      x_ann <- 10^xy_lim[1]; y_ann <- 10^xy_lim[2]
    } else {
      x_ann <- xy_lim[1]; y_ann <- xy_lim[2]
    }
  } else {
    x_ann <- min(d$x, na.rm = TRUE); y_ann <- max(d$y, na.rm = TRUE)
  }
  p <- ggplot(d, aes(x = x, y = y)) +
    geom_point(alpha = 0.15, size = 0.5, color = RICE_BLUE) +
    geom_abline(slope = 1, intercept = 0, color = HHD_ORANGE,
                linetype = "dashed", linewidth = 0.5) +
    annotate("label", x = x_ann, y = y_ann, label = ann,
             hjust = 0, vjust = 1, size = 2.6,
             fill = "white", alpha = 0.85,
             color = RICE_BLUE_DARK, family = "mono") +
    labs(
      title    = NULL,
      subtitle = NULL,
      x = x_lab, y = y_lab
    ) +
    theme_minimal(base_size = 9) +
    theme(plot.title    = element_text(face = "bold",
                                        color = RICE_BLUE_DARK),
          plot.subtitle = element_text(color = NEUTRAL_GRAY, size = 8))
  if (log_log) {
    p <- p + scale_x_log10() + scale_y_log10()
    if (!is.null(xy_lim))
      p <- p + coord_fixed(xlim = 10^c(xy_lim[1], xy_lim[2]),
                            ylim = 10^c(xy_lim[1], xy_lim[2]))
    else
      p <- p + coord_fixed()
  } else {
    if (!is.null(xy_lim))
      p <- p + coord_fixed(xlim = c(xy_lim[1], xy_lim[2]),
                            ylim = c(xy_lim[1], xy_lim[2]))
    else
      p <- p + coord_fixed()
  }
  p
}

# Rt scatter matrix (todos los pares — para 8 fits son 28 paneles)
cat(sprintf("  Building fig9_Rt_scatter_all_pairs.pdf (%d pares) ...\n",
            length(combos)))
n_cols_scatter <- ifelse(length(combos) > 15, 4, 3)
n_rows_scatter <- ceiling(length(combos) / n_cols_scatter)
scatter_height <- max(8, 2.4 * n_rows_scatter)

rt_all <- unlist(lapply(combos, function(cb) {
  c(rt_mean[[cb[1]]], rt_mean[[cb[2]]])
}))
rt_xy_lim <- c(0, ceiling(max(rt_all, na.rm = TRUE) * 10) / 10)
rt_panels <- lapply(combos, function(cb) {
  make_scatter_panel(rt_mean, cb[1], cb[2], log_log = FALSE,
                     title_extra = "R_t posterior median",
                     xy_lim = rt_xy_lim)
})
rt_combo_grob <- gridExtra::arrangeGrob(
  grobs = rt_panels, ncol = n_cols_scatter)
ggsave(file.path(OUT_DIR, "fig9_Rt_scatter_all_pairs_AR2.pdf"),
       rt_combo_grob, width = 14, height = scatter_height,
       limitsize = FALSE)
cat("  -> fig9_Rt_scatter_all_pairs.pdf\n")

# It scatter matrix (log-log)
cat(sprintf("  Building fig10_It_scatter_all_pairs.pdf (%d pares) ...\n",
            length(combos)))
it_all <- unlist(lapply(combos, function(cb) {
  c(it_mean[[cb[1]]], it_mean[[cb[2]]])
}))
it_all <- it_all[is.finite(it_all) & it_all > 0]
it_xy_lim <- c(floor(log10(min(it_all)) * 10) / 10,
                ceiling(log10(max(it_all)) * 10) / 10)
it_panels <- lapply(combos, function(cb) {
  make_scatter_panel(it_mean, cb[1], cb[2], log_log = TRUE,
                     title_extra = "I_t posterior median, log10",
                     xy_lim = it_xy_lim)
})
it_combo_grob <- gridExtra::arrangeGrob(
  grobs = it_panels, ncol = n_cols_scatter)
ggsave(file.path(OUT_DIR, "fig10_It_scatter_all_pairs_AR2.pdf"),
       it_combo_grob, width = 14, height = scatter_height,
       limitsize = FALSE)
cat("  -> fig10_It_scatter_all_pairs.pdf\n\n")

# ============================================================================
# 10 . Diagnosticos por bloque (Rhat y ESS) — estilo paper Fig 6/7 del .qmd
# ============================================================================
cat("== 10 . Diagnosticos por bloque (Rhat y ESS) ==\n")

# Construir long df con rhat y ess_bulk por variable, taggeado por block
rhat_long <- imap_dfr(fits, function(fit, lab) {
  s <- fit$summary
  s$block <- case_when(
    grepl("^Rt\\[",           s$variable) ~ "Rt",
    grepl("^It\\[",           s$variable) ~ "It",
    grepl("^beta\\[",         s$variable) ~ "beta",
    grepl("^log_beta\\[",     s$variable) ~ "log_beta",
    grepl("^sigma_obs_pl\\[", s$variable) ~ "sigma_obs",
    s$variable %in% c("rho")              ~ "rho",
    TRUE                                  ~ "other"
  )
  s$fit <- lab
  s |> filter(block %in% c("Rt", "It", "beta", "rho"))
})

# Tabla resumen Rhat y ESS por (fit, block)
rhat_summary <- rhat_long |>
  group_by(fit, block) |>
  summarise(
    n              = n(),
    rhat_max       = round(max(rhat,        na.rm = TRUE), 4),
    rhat_med       = round(median(rhat,     na.rm = TRUE), 4),
    pct_above_1.05 = round(mean(rhat > 1.05, na.rm = TRUE) * 100, 2),
    ess_min        = round(min(ess_bulk,    na.rm = TRUE)),
    ess_med        = round(median(ess_bulk, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  arrange(block, fit)
print(rhat_summary, row.names = FALSE)
write.csv(rhat_summary, file.path(OUT_DIR, "rhat_ess_summary_by_block_AR2.csv"),
          row.names = FALSE)
cat("  -> rhat_ess_summary_by_block.csv\n")

# Rhat violin/box plot por bloque (estilo paper 06_rhat_by_block.pdf)
# El block "rho" tiene solo 1 fila por fit rho-est (no se puede violin).
# Lo dejamos fuera de las violins; aparece en table tableS_phi_posterior.csv.

# ESS bulk violin/box plot por bloque
p_ess <- ggplot(rhat_long |> filter(block != "rho"),
                  aes(x = fit, y = ess_bulk, fill = fit)) +
  geom_hline(yintercept = 400, linetype = "dashed",
             color = HHD_ORANGE_DARK, linewidth = 0.5) +
  geom_violin(alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.18, outlier.size = 0.4, fill = "white") +
  scale_fill_manual(values = COLORS_4, guide = "none") +
  scale_x_discrete(labels = fit_label_parser) +
  scale_y_log10(labels = comma) +
  facet_wrap(~ block, scales = "free_y", nrow = 1) +
  labs(
    title    = NULL,
    subtitle = NULL,
    x = NULL, y = "ESS bulk (log10)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", color = RICE_BLUE_DARK),
        plot.subtitle = element_text(color = NEUTRAL_GRAY, size = 9),
        strip.text    = element_text(face = "bold"),
        axis.text.x   = element_text(angle = 35, hjust = 1, size = 9))

ggsave(file.path(OUT_DIR, "fig12_ess_by_block_AR2.pdf"),
       p_ess, width = 13, height = 5)
cat("  -> fig12_ess_by_block.pdf\n")

# Lollipop estilo paper: Rhat worst-case per plant, by fit y bloque
# (mas compacto que el paper original — sin sweep de rho)
rhat_lolli <- rhat_long |>
  filter(block %in% c("Rt", "It", "beta")) |>
  mutate(p = case_when(
    block == "beta" ~ as.integer(str_match(variable, "beta\\[(\\d+)\\]")[, 2]),
    TRUE            ~ as.integer(str_match(variable, "\\[(\\d+),\\d+\\]")[, 2])
  )) |>
  filter(!is.na(p), p >= 1, p <= P) |>
  group_by(fit, block, p) |>
  summarise(rhat_worst = max(rhat, na.rm = TRUE), .groups = "drop") |>
  mutate(
    plant     = factor(plant_names[p],
                       levels = plant_names,
                       labels = plant_labels[plant_names]),
    converged = rhat_worst <= 1.05
  )

y_max <- max(rhat_lolli$rhat_worst, na.rm = TRUE)

p_lolli <- ggplot(rhat_lolli, aes(x = plant, y = rhat_worst, color = fit)) +
  geom_segment(aes(xend = plant, yend = 1),
               linewidth = 0.35, alpha = 0.6) +
  geom_point(aes(fill = converged), size = 2.2, shape = 21, stroke = 0.6) +
  geom_hline(yintercept = 1.05, linetype = "dashed",
             color = HHD_ORANGE_DARK, linewidth = 0.4) +
  geom_hline(yintercept = 1.01, linetype = "dashed",
             color = "darkgreen", linewidth = 0.4) +
  geom_hline(yintercept = 1.00, color = "grey85", linewidth = 0.3) +
  scale_color_manual(values = COLORS_4, name = NULL,
                     labels = fit_label_parser) +
  scale_fill_manual(
    values = c(`TRUE` = "white", `FALSE` = HHD_ORANGE_DARK),
    labels = c(`TRUE` = "Converged", `FALSE` = "Not converged"),
    name   = expression(hat(R) <= 1.05)
  ) +
  facet_grid(fit ~ block,
             labeller = labeller(
               fit   = as_labeller(title_label_parse, label_parsed),
               block = block_labeller
             )) +
  coord_cartesian(ylim = c(0.99, y_max + 0.005)) +
  labs(
    title    = NULL,
    subtitle = NULL,
    x = NULL, y = expression(hat(R) * " (worst per plant)")
  ) +
  theme_minimal(base_size = 9) +
  theme(plot.title       = element_text(face = "bold",
                                          color = RICE_BLUE_DARK),
        plot.subtitle    = element_text(color = NEUTRAL_GRAY, size = 8),
        axis.text.x      = element_text(angle = 60, hjust = 1, size = 7),
        strip.text       = element_text(face = "bold", size = 8),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom")

ggsave(file.path(OUT_DIR, "fig13_rhat_lollipop_AR2.pdf"),
       p_lolli, width = 13, height = max(6, length(fits) * 1.8))
cat("  -> fig13_rhat_lollipop.pdf\n\n")

# ============================================================================
# 11 . CRPS — Continuous Ranked Probability Score
# ----------------------------------------------------------------------------
# Para cada fit, CRPS contra las observaciones y_raw (NO los y_filtrados).
# El modelo observacional es:
#   y_it ~ Lognormal(log(psi_it + 1e-6) - 0.5 * sigma_it^2, sigma_it)
# donde psi_it = beta_i * sum_{tau=0}^{D-1} d_tau * It[i, t-tau].
#
# Para cada (planta, semana) observada (y_raw no NA), calculamos:
#   - posterior predictive samples de y_it
#   - CRPS via scoringRules::crps_sample()
# Promediar -> 1 numero por fit. Lower CRPS = mejor ajuste predictivo.
#
# Ademas CRPS en log10 para tratar plantas con igual peso (las grandes
# dominan en escala original).
# ============================================================================
cat("== 11 . CRPS (posterior predictive vs y_raw) ==\n")

# ---- helpers --------------------------------------------------------------
get_It_array <- function(fit, draws_idx) {
  # Devuelve array [length(draws_idx), P, T] con It[d, p, t]
  it_mat <- get_draws_matrix(fit, "It")[draws_idx, , drop = FALSE]
  # Stan: It[p, t] -> columna p + (t-1)*P en el draws_matrix
  T_loc <- ncol(it_mat) / P
  array(it_mat, dim = c(length(draws_idx), P, T_loc))
}

compute_psi_samples <- function(fit, n_thin = 500) {
  beta_mat <- get_draws_matrix(fit, "beta")  # [n_draws, P]
  n_draws  <- nrow(beta_mat)
  draws_idx <- round(seq(1, n_draws, length.out = min(n_thin, n_draws)))
  beta_sub <- beta_mat[draws_idx, , drop = FALSE]      # [n_thin, P]
  It_arr   <- get_It_array(fit, draws_idx)             # [n_thin, P, T]
  T_loc <- dim(It_arr)[3]

  d_vec <- as.numeric(fit$d_vec)
  D <- length(d_vec)

  # Convolucion vectorizada: psi[d, p, t] = beta[d, p] * sum_tau d_tau It[d, p, t-tau]
  psi_arr <- array(NA_real_,
                    dim = c(length(draws_idx), P, T_loc))
  for (t in seq_len(T_loc)) {
    # tau valido: 0..D-1, con t-tau >= 1
    max_tau <- min(D - 1, t - 1)
    if (max_tau < 0) next
    conv <- matrix(0, nrow = length(draws_idx), ncol = P)
    for (tau in 0:max_tau) {
      conv <- conv + d_vec[tau + 1] * It_arr[, , t - tau]
    }
    psi_arr[, , t] <- beta_sub * conv
  }
  list(psi = psi_arr, n_thin = length(draws_idx), draws_idx = draws_idx)
}

get_sigma_obs_arr <- function(fit, draws_idx) {
  # Determine T from It summary
  it_summary <- fit$summary[grepl("^It\\[", fit$summary$variable), ]
  T_loc <- length(it_summary$mean) / P

  # 1. SSM con sigma_filt guardado -> array [n_thin, P, T]
  if (!is.null(fit$sigma_filt) && is.matrix(fit$sigma_filt) &&
      all(dim(fit$sigma_filt) == c(P, T_loc))) {
    sf <- fit$sigma_filt
    arr <- array(NA_real_, dim = c(length(draws_idx), P, T_loc))
    for (d in seq_along(draws_idx)) arr[d, , ] <- sf
    return(arr)
  }

  # 2. Joint con draws de sigma_obs_pl -> usar draws
  sm <- tryCatch(get_draws_matrix(fit, "sigma_obs_pl"),
                  error = function(e) NULL)
  if (!is.null(sm) && ncol(sm) == P) {
    sm <- sm[draws_idx, , drop = FALSE]
    arr <- array(NA_real_, dim = c(length(draws_idx), P, T_loc))
    for (t in seq_len(T_loc)) arr[, , t] <- sm
    return(arr)
  }

  # 3. Fallback: posterior means de sigma_obs_pl en summary
  sigma_summary <- fit$summary[grepl("^sigma_obs_pl\\[", fit$summary$variable), ]
  if (nrow(sigma_summary) == P) {
    sm_means <- matrix(sigma_summary$mean,
                       nrow = length(draws_idx), ncol = P, byrow = TRUE)
    arr <- array(NA_real_, dim = c(length(draws_idx), P, T_loc))
    for (t in seq_len(T_loc)) arr[, , t] <- sm_means
    return(arr)
  }

  # 4. No hay sigma de ningun lado -> NULL para indicar skip
  warning(sprintf("No sigma_obs encontrado en el fit (ni sigma_filt, ni draws, ni summary). Skip CRPS."))
  NULL
}

compute_predictive_samples <- function(fit, n_thin = 500) {
  ps <- compute_psi_samples(fit, n_thin)
  psi_arr <- ps$psi
  sig_arr <- get_sigma_obs_arr(fit, ps$draws_idx)
  if (is.null(sig_arr)) return(NULL)   # no sigma -> skip
  if (!all(dim(psi_arr) == dim(sig_arr))) {
    warning("Dimensiones de psi y sigma no coinciden, skip")
    return(NULL)
  }
  # Sample y_it from Lognormal(log(psi+eps) - 0.5*sig^2, sig)
  y_pred <- array(NA_real_, dim = dim(psi_arr))
  finite_idx <- is.finite(psi_arr) & psi_arr > 0 & is.finite(sig_arr) & sig_arr > 0
  mu_log <- log(psi_arr[finite_idx] + 1e-6) - 0.5 * sig_arr[finite_idx]^2
  y_pred[finite_idx] <- exp(stats::rnorm(sum(finite_idx),
                                          mean = mu_log,
                                          sd   = sig_arr[finite_idx]))
  y_pred
}

compute_crps_matrices <- function(y_pred, y_obs) {
  # y_pred: [n_thin, P, T_pred], y_obs: [P, T_obs] (idealmente matriz P x T)
  # Validacion defensiva (algunos fits viejos guardaron y_raw con formato distinto)
  if (is.null(y_obs) || !is.matrix(y_obs)) {
    warning("y_obs no es matriz P x T; skip CRPS para este fit")
    return(NULL)
  }
  if (nrow(y_obs) != P) {
    warning(sprintf("y_obs tiene %d filas, esperaba P=%d; skip CRPS",
                     nrow(y_obs), P))
    return(NULL)
  }
  T_use <- min(dim(y_pred)[3], ncol(y_obs))
  crps_orig <- matrix(NA_real_, P, T_use)
  crps_log  <- matrix(NA_real_, P, T_use)
  for (t in seq_len(T_use)) {
    for (p in seq_len(P)) {
      yo <- y_obs[p, t]
      # Defensa: yo debe ser un escalar finito positivo
      if (length(yo) != 1 || !is.finite(yo) || yo <= 0) next
      yp <- y_pred[, p, t]
      yp <- yp[is.finite(yp) & yp > 0]
      if (length(yp) < 20) next
      crps_orig[p, t] <- scoringRules::crps_sample(y = yo, dat = yp)
      crps_log[p, t]  <- scoringRules::crps_sample(y = log10(yo),
                                                    dat = log10(yp))
    }
  }
  list(orig = crps_orig, log10 = crps_log)
}

# ---- Enrichment: presta y_raw, d_vec, sigma_filt de fits que los tienen ----
# Los fits viejos del rho-sweep (A_ssm_fix, A_direct_fix) no guardaron y_raw
# ni d_vec (son del dataset, no del modelo). Los reconstruimos prestandolos.
ref_yraw_dvec <- NULL
for (lab in names(fits)) {
  if (!is.null(fits[[lab]]$y_raw) && is.matrix(fits[[lab]]$y_raw) &&
      !is.null(fits[[lab]]$d_vec)) {
    ref_yraw_dvec <- fits[[lab]]
    cat(sprintf("  Referencia y_raw+d_vec: %s\n", lab))
    break
  }
}
ref_sigma_filt <- NULL
for (lab in names(fits)) {
  if (!is.null(fits[[lab]]$sigma_filt) && is.matrix(fits[[lab]]$sigma_filt)) {
    ref_sigma_filt <- fits[[lab]]$sigma_filt
    cat(sprintf("  Referencia sigma_filt: %s\n", lab))
    break
  }
}

for (lab in names(fits)) {
  if (is.null(fits[[lab]]$y_raw) && !is.null(ref_yraw_dvec)) {
    fits[[lab]]$y_raw <- ref_yraw_dvec$y_raw
    cat(sprintf("  -> %s: inyectado y_raw del referencia\n", lab))
  }
  if (is.null(fits[[lab]]$d_vec) && !is.null(ref_yraw_dvec)) {
    fits[[lab]]$d_vec <- ref_yraw_dvec$d_vec
  }
  # Para fits SSM (los _fix con use_filter=1) que no guardaron sigma_filt
  if (grepl("ssm", lab) && is.null(fits[[lab]]$sigma_filt) &&
      !is.null(ref_sigma_filt)) {
    fits[[lab]]$sigma_filt <- ref_sigma_filt
    cat(sprintf("  -> %s: inyectado sigma_filt del referencia (SSM)\n", lab))
  }
}
cat("\n")

# ---- correr CRPS por cada fit --------------------------------------------
N_THIN <- 500   # subset de draws para CRPS (precision adecuada, costo razonable)
crps_results <- list()
for (lab in names(fits)) {
  cat(sprintf("  Computing CRPS for %-18s ", lab))
  t0 <- Sys.time()
  fit <- fits[[lab]]
  y_pred <- tryCatch(
    compute_predictive_samples(fit, n_thin = N_THIN),
    error = function(e) {
      cat(sprintf("FAIL: %s\n", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(y_pred)) {
    cat("(skipped: no predictive samples)\n")
    next
  }
  cm <- tryCatch(
    compute_crps_matrices(y_pred, fit$y_raw),
    error = function(e) {
      cat(sprintf("FAIL en compute_crps_matrices: %s\n",
                   conditionMessage(e)))
      NULL
    }
  )
  if (is.null(cm)) {
    cat("(skipped: y_raw malformado o sin valores validos)\n")
    next
  }
  crps_results[[lab]] <- cm
  cat(sprintf("(%.1f sec)\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

if (length(crps_results) == 0) {
  cat("  WARN: ningun fit pudo calcular CRPS — saltando seccion\n\n")
}

if (length(crps_results) > 0) {

# ---- agregaciones ---------------------------------------------------------
# Por fit: mean CRPS overall (ambas escalas) + mediana
crps_overall <- imap_dfr(crps_results, function(cm, lab) {
  data.frame(
    fit            = lab,
    fit_label      = title_label_plain[[lab]],
    n_obs          = sum(!is.na(cm$orig)),
    crps_mean      = mean(cm$orig,  na.rm = TRUE),
    crps_median    = median(cm$orig, na.rm = TRUE),
    crps_log_mean  = mean(cm$log10, na.rm = TRUE),
    crps_log_median= median(cm$log10, na.rm = TRUE)
  )
})
print(crps_overall, row.names = FALSE, digits = 4)
write.csv(crps_overall, file.path(OUT_DIR, "crps_summary_AR2.csv"),
          row.names = FALSE)
cat("  -> crps_summary.csv\n")

# Por (fit, planta): breakdown per-plant
crps_by_plant <- imap_dfr(crps_results, function(cm, lab) {
  data.frame(
    fit            = lab,
    plant          = plant_names,
    plant_lab      = plant_labels,
    n_obs          = rowSums(!is.na(cm$orig)),
    crps_mean      = rowMeans(cm$orig,  na.rm = TRUE),
    crps_log_mean  = rowMeans(cm$log10, na.rm = TRUE)
  )
})
write.csv(crps_by_plant, file.path(OUT_DIR, "crps_by_plant_AR2.csv"),
          row.names = FALSE)
cat("  -> crps_by_plant.csv\n")

# ---- Figura 14: CRPS overall (orig + log10) -------------------------------
crps_overall_long <- crps_overall |>
  select(fit, crps_mean, crps_log_mean) |>
  pivot_longer(cols = c(crps_mean, crps_log_mean),
               names_to = "scale", values_to = "crps") |>
  mutate(scale_lab = ifelse(scale == "crps_mean",
                             "CRPS (B gc/day)",
                             "CRPS (log10 scale)"))

p_fig14 <- ggplot(crps_overall_long,
                  aes(x = reorder(fit, crps), y = crps, fill = fit)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.3g", crps)),
            hjust = -0.1, size = 3) +
  scale_fill_manual(values = COLORS_4, guide = "none") +
  scale_x_discrete(labels = fit_label_parser) +
  coord_flip() +
  facet_wrap(~ scale_lab, scales = "free_x", ncol = 2) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL, y = "Mean CRPS"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", color = RICE_BLUE_DARK),
        plot.subtitle = element_text(color = NEUTRAL_GRAY),
        strip.text    = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "fig14_crps_summary_AR2.pdf"),
       p_fig14, width = 11, height = 5)
cat("  -> fig14_crps_summary.pdf\n")


# ---- Figura 14: CRPS overall (orig + log10) -------------------------------

# ---- Figura 15: CRPS por planta -------------------------------------------
crps_by_plant_plot <- crps_by_plant |>
  mutate(plant = factor(plant_lab,
                         levels = plant_labels[plant_names]))

p_fig15 <- ggplot(crps_by_plant_plot,
                   aes(x = plant, y = crps_log_mean,
                       fill = fit, color = fit)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7, alpha = 0.85) +
  scale_fill_manual(values = COLORS_4, name = NULL,
                    labels = fit_label_parser) +
  scale_color_manual(values = COLORS_4, guide = "none") +
  labs(
    title    = NULL,
    subtitle = NULL,
    x = NULL, y = expression("CRPS, log"[10]*" scale")
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title    = element_text(face = "bold", color = RICE_BLUE_DARK),
        plot.subtitle = element_text(color = NEUTRAL_GRAY),
        axis.text.x   = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(OUT_DIR, "fig15_crps_by_plant_AR2.pdf"),
       p_fig15, width = 14, height = 6)
cat("  -> fig15_crps_by_plant.pdf\n\n")

# ============================================================================
# 11b . SANITY CHECK CRPS — solo SSM fits, evaluados vs y_filt (no y_raw)
# ----------------------------------------------------------------------------
# Razon: SSM se entrena contra y_filt + sigma_filt (Stage-1 Kalman). Su
# posterior predictive es Lognormal(log psi - 0.5*sigma_filt^2, sigma_filt).
# Si SSM es internamente consistente, deberia tener CRPS MUY bajo vs y_filt
# (mucho mejor que vs y_raw). Eso confirma que la diferencia que vimos contra
# Joint en CRPS principal viene del mismatch target (y_filt entrenamiento vs
# y_raw evaluacion), no de un problema fundamental con SSM.
# ============================================================================
cat("== 11b . SANITY CHECK: SSM fits vs y_filt (target propio del SSM) ==\n")

# Necesitamos y_filt — lo buscamos en cualquier fit que lo tenga
ref_yfilt <- NULL
for (lab in names(fits)) {
  if (!is.null(fits[[lab]]$y_filt) && is.matrix(fits[[lab]]$y_filt)) {
    ref_yfilt <- fits[[lab]]$y_filt
    cat(sprintf("  Usando y_filt de %s como target\n", lab))
    break
  }
}

if (is.null(ref_yfilt)) {
  cat("  WARN: ningun fit tiene y_filt — sanity check imposible, skip\n\n")
} else {

  # Solo SSM fits (Joint no aplica — su target es y_raw)
  ssm_labs <- grep("^A_ssm_", names(fits), value = TRUE)
  cat(sprintf("  SSM fits a evaluar: %s\n",
              paste(ssm_labs, collapse=", ")))

  crps_sanity <- list()
  for (lab in ssm_labs) {
    cat(sprintf("  Computing sanity CRPS for %-18s ", lab))
    t0 <- Sys.time()
    fit <- fits[[lab]]
    y_pred <- tryCatch(
      compute_predictive_samples(fit, n_thin = N_THIN),
      error = function(e) {
        cat(sprintf("FAIL: %s\n", conditionMessage(e)))
        NULL
      }
    )
    if (is.null(y_pred)) {
      cat("(skipped)\n"); next
    }
    cm <- tryCatch(
      compute_crps_matrices(y_pred, ref_yfilt),   # <- target = y_filt
      error = function(e) { cat(sprintf("FAIL: %s\n",
                                       conditionMessage(e))); NULL }
    )
    if (is.null(cm)) {
      cat("(skipped)\n"); next
    }
    crps_sanity[[lab]] <- cm
    cat(sprintf("(%.1f sec)\n",
                as.numeric(difftime(Sys.time(), t0, units="secs"))))
  }

  if (length(crps_sanity) > 0) {
    # Resumen lado a lado: CRPS vs y_raw (principal) y CRPS vs y_filt (sanity)
    sanity_overall <- imap_dfr(crps_sanity, function(cm, lab) {
      data.frame(
        fit             = lab,
        fit_label       = title_label_plain[[lab]],
        n_obs           = sum(!is.na(cm$orig)),
        crps_mean_sanity      = mean(cm$orig,  na.rm = TRUE),
        crps_log_mean_sanity  = mean(cm$log10, na.rm = TRUE)
      )
    })
    # Merge con el CRPS principal para mostrar comparacion
    main_for_ssm <- crps_overall |>
      filter(fit %in% ssm_labs) |>
      select(fit, crps_mean_main = crps_mean,
             crps_log_mean_main = crps_log_mean)
    sanity_compare <- sanity_overall |>
      left_join(main_for_ssm, by = "fit") |>
      mutate(
        ratio_orig = crps_mean_main / crps_mean_sanity,
        ratio_log  = crps_log_mean_main / crps_log_mean_sanity
      ) |>
      select(fit, fit_label,
             crps_mean_main, crps_mean_sanity, ratio_orig,
             crps_log_mean_main, crps_log_mean_sanity, ratio_log)
    print(sanity_compare, row.names = FALSE, digits = 4)
    write.csv(sanity_compare,
              file.path(OUT_DIR, "crps_ssm_sanity_check_AR2.csv"),
              row.names = FALSE)
    cat("  -> crps_ssm_sanity_check.csv\n\n")
    cat("INTERPRETACION:\n")
    cat("  ratio_orig (y_raw/y_filt) >> 1 -> SSM ajusta bien y_filt pero la\n")
    cat("  inflacion al cambiar target a y_raw confirma el mismatch.\n\n")
  } else {
    cat("  Ningun fit SSM pudo evaluar sanity. Skip.\n\n")
  }
}

}  # cierre del if (length(crps_results) > 0)

cat("============================================================\n")
cat("  DONE — outputs en:\n")
cat(sprintf("  %s\n", OUT_DIR))
cat("============================================================\n")
cat("\nArchivos generados:\n")
cat("  Figuras estilo paper:\n")
cat("    fig6_phi_posterior.pdf            (Figure 6 paper)\n")
cat("    fig7_Rt_by_plant.pdf              (Figure 7 paper, 4 paneles)\n")
cat("    fig8_It_by_plant.pdf              (Figure 8 paper, 4 paneles)\n")
cat("  Correlaciones entre todos los fits:\n")
cat("    fig9_Rt_scatter_all_pairs.pdf     (6 pares, scatter Rt)\n")
cat("    fig10_It_scatter_all_pairs.pdf    (6 pares, scatter It log-log)\n")
cat("  Diagnosticos:\n")
cat("    fig12_ess_by_block.pdf            (ESS bulk violin/box por bloque)\n")
cat("    fig13_rhat_lollipop.pdf           (Rhat worst-case per plant)\n")
cat("  Tablas:\n")
cat("    table1_beta_posterior.csv         (Tabla 1 paper)\n")
cat("    tableS_phi_posterior.csv          (phi posterior)\n")
cat("    agreement_all_pairs.csv           (cor R, cor I log10, ratios)\n")
cat("    rhat_ess_summary_by_block.csv     (resumen Rhat/ESS por bloque)\n")
cat("    convergence_diagnostic.csv\n")
cat("    loo_compare.csv\n")
cat("  CSVs largos:\n")
cat("    Rt_posterior_all.csv\n")
cat("    It_posterior_all.csv\n")
cat("    beta_posterior_long.csv\n")
cat("  Scoring:\n")
cat("    crps_summary.csv                  (mean CRPS per fit)\n")
cat("    crps_by_plant.csv                 (per-plant breakdown)\n")
cat("    fig15_crps_by_plant.pdf           (CRPS by plant + fit)\n\n")
