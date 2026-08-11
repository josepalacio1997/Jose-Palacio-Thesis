#!/usr/bin/env bash
# =============================================================================
# fix_scatter_panel_sizes.sh
# -----------------------------------------------------------------------------
# Makes all 6 scatter panels in fig9/fig10 the same physical size by computing
# global x/y limits across all data and applying them uniformly.
#
# Without this, coord_equal() lets each panel resize to fit its own data range,
# so panels with larger values (e.g., Data; rho-hat reaching R_t > 3) end up
# physically larger than those reaching only R_t ~ 2.
#
# Run from the repo root:
#   bash "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/fix_scatter_panel_sizes.sh"
# =============================================================================

set -e

cd "$HOME/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW"

FILES=(
  "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/analyze_models_A_AR1.R"
  "Modelos AR con rho est. y rho=0/paper_pipeline/08_analysis_AR2/analyze_models_A_AR2.R"
)

DATE=$(date +%Y%m%d_panelfix)
for f in "${FILES[@]}"; do
  cp "$f" "$f.bak.$DATE"
done

python3 << 'PYEOF'
from pathlib import Path

files = [
    "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/analyze_models_A_AR1.R",
    "Modelos AR con rho est. y rho=0/paper_pipeline/08_analysis_AR2/analyze_models_A_AR2.R",
]

# Replace make_scatter_panel: add common-limit argument, use coord_fixed with
# explicit limits so all panels render at identical size.
old_block = '''make_scatter_panel <- function(df, a, b, log_log = FALSE, title_extra = "") {
  d <- data.frame(x = df[[a]], y = df[[b]]) |>
    filter(is.finite(x) & is.finite(y))
  if (log_log) d <- filter(d, x > 0, y > 0)
  if (log_log) {
    xx <- log10(d$x); yy <- log10(d$y)
  } else {
    xx <- d$x; yy <- d$y
  }
  r    <- cor(xx, yy)
  rmse <- sqrt(mean((xx - yy)^2))
  ann  <- sprintf("r = %.3f\\nRMSE = %.3f", r, rmse)
  # Map fit IDs to formatted axis labels
  x_lab <- title_label[[a]]; if (is.null(x_lab)) x_lab <- a
  y_lab <- title_label[[b]]; if (is.null(y_lab)) y_lab <- b
  # Top-left annotation in DATA coords
  x_ann <- min(d$x, na.rm = TRUE)
  y_ann <- max(d$y, na.rm = TRUE)
  p <- ggplot(d, aes(x = x, y = y)) +
    geom_point(alpha = 0.15, size = 0.5, color = RICE_BLUE) +
    geom_abline(slope = 1, intercept = 0, color = HHD_ORANGE,
                linetype = "dashed", linewidth = 0.5) +
    annotate("label", x = x_ann, y = y_ann, label = ann,
             hjust = 0, vjust = 1, size = 2.6,
             label.size = 0, fill = "white", alpha = 0.85,
             color = RICE_BLUE_DARK, family = "mono") +
    coord_equal() +
    labs(x = x_lab, y = y_lab) +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(face = "bold", color = RICE_BLUE_DARK),
          plot.subtitle = element_text(color = NEUTRAL_GRAY, size = 8))
  if (log_log) p <- p + scale_x_log10() + scale_y_log10()
  p
}'''

new_block = '''make_scatter_panel <- function(df, a, b, log_log = FALSE, title_extra = "",
                                xy_lim = NULL) {
  d <- data.frame(x = df[[a]], y = df[[b]]) |>
    filter(is.finite(x) & is.finite(y))
  if (log_log) d <- filter(d, x > 0, y > 0)
  if (log_log) {
    xx <- log10(d$x); yy <- log10(d$y)
  } else {
    xx <- d$x; yy <- d$y
  }
  r    <- cor(xx, yy)
  rmse <- sqrt(mean((xx - yy)^2))
  ann  <- sprintf("r = %.3f\\nRMSE = %.3f", r, rmse)
  x_lab <- title_label[[a]]; if (is.null(x_lab)) x_lab <- a
  y_lab <- title_label[[b]]; if (is.null(y_lab)) y_lab <- b
  # Annotation placed in the upper-left corner using the shared limits
  if (!is.null(xy_lim)) {
    x_ann <- xy_lim[1]; y_ann <- xy_lim[2]
  } else {
    x_ann <- min(d$x, na.rm = TRUE); y_ann <- max(d$y, na.rm = TRUE)
  }
  p <- ggplot(d, aes(x = x, y = y)) +
    geom_point(alpha = 0.15, size = 0.5, color = RICE_BLUE) +
    geom_abline(slope = 1, intercept = 0, color = HHD_ORANGE,
                linetype = "dashed", linewidth = 0.5) +
    annotate("label", x = x_ann, y = y_ann, label = ann,
             hjust = 0, vjust = 1, size = 2.6,
             label.size = 0, fill = "white", alpha = 0.85,
             color = RICE_BLUE_DARK, family = "mono") +
    labs(x = x_lab, y = y_lab) +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(face = "bold", color = RICE_BLUE_DARK),
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
}'''

# rt_panels / it_panels call sites: pass the common limits
old_rt_call = '''rt_panels <- lapply(combos, function(cb) {
  make_scatter_panel(rt_mean, cb[1], cb[2], log_log = FALSE,
                     title_extra = "R_t posterior median")
})'''
new_rt_call = '''rt_all <- unlist(lapply(combos, function(cb) {
  c(rt_mean[[cb[1]]], rt_mean[[cb[2]]])
}))
rt_xy_lim <- c(0, ceiling(max(rt_all, na.rm = TRUE) * 10) / 10)
rt_panels <- lapply(combos, function(cb) {
  make_scatter_panel(rt_mean, cb[1], cb[2], log_log = FALSE,
                     title_extra = "R_t posterior median",
                     xy_lim = rt_xy_lim)
})'''

old_it_call = '''it_panels <- lapply(combos, function(cb) {
  make_scatter_panel(it_mean, cb[1], cb[2], log_log = TRUE,
                     title_extra = "I_t posterior median, log10")
})'''
new_it_call = '''it_all <- unlist(lapply(combos, function(cb) {
  c(it_mean[[cb[1]]], it_mean[[cb[2]]])
}))
it_all <- it_all[is.finite(it_all) & it_all > 0]
it_xy_lim <- c(floor(log10(min(it_all)) * 10) / 10,
                ceiling(log10(max(it_all)) * 10) / 10)
it_panels <- lapply(combos, function(cb) {
  make_scatter_panel(it_mean, cb[1], cb[2], log_log = TRUE,
                     title_extra = "I_t posterior median, log10",
                     xy_lim = it_xy_lim)
})'''

for f in files:
    p = Path(f)
    s = p.read_text()
    if old_block in s:
        s = s.replace(old_block, new_block)
        print(f"  patched make_scatter_panel in {f.split('/')[-1]}")
    else:
        print(f"  WARN: make_scatter_panel pattern not matched in {f.split('/')[-1]}")
    if old_rt_call in s:
        s = s.replace(old_rt_call, new_rt_call)
        print(f"    + rt_panels call site")
    if old_it_call in s:
        s = s.replace(old_it_call, new_it_call)
        print(f"    + it_panels call site")
    p.write_text(s)

print()
print("Done. Now re-run the analyses to regenerate the figures:")
print()
print("  Rscript \"Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/analyze_models_A_AR1.R\"")
print("  Rscript \"Modelos AR con rho est. y rho=0/paper_pipeline/08_analysis_AR2/analyze_models_A_AR2.R\"")
PYEOF
