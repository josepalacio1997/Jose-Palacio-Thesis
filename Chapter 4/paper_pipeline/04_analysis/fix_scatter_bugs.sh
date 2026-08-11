#!/usr/bin/env bash
# =============================================================================
# fix_scatter_bugs.sh
# -----------------------------------------------------------------------------
# Fixes three bugs in the scatter plot generation:
#   1. AR2 make_scatter_panel doesn't accept xy_lim → add it
#   2. fig10 (log-log) annotation invisible because x_ann was set to the log10
#      value but scale_x_log10 expects natural scale → set 10^xy_lim
#   3. annotate(..., label.size = 0) warning on older ggplot2 → drop param
# =============================================================================

set -e
cd "$HOME/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW"

FILES=(
  "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/analyze_models_A_AR1.R"
  "Modelos AR con rho est. y rho=0/paper_pipeline/08_analysis_AR2/analyze_models_A_AR2.R"
)

DATE=$(date +%Y%m%d_scatterbugs)
for f in "${FILES[@]}"; do
  cp "$f" "$f.bak.$DATE"
done

python3 << 'PYEOF'
from pathlib import Path

files = [
    "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/analyze_models_A_AR1.R",
    "Modelos AR con rho est. y rho=0/paper_pipeline/08_analysis_AR2/analyze_models_A_AR2.R",
]

# Common new annotate (no label.size, correct log-log coords)
NEW_ANNOTATE = '''  # Top-left annotation: data-scale coords for log_log, raw for linear
  if (!is.null(xy_lim)) {
    if (log_log) {
      x_ann <- 10^xy_lim[1]; y_ann <- 10^xy_lim[2]
    } else {
      x_ann <- xy_lim[1]; y_ann <- xy_lim[2]
    }
  } else {
    x_ann <- min(d$x, na.rm = TRUE); y_ann <- max(d$y, na.rm = TRUE)
  }'''

NEW_ANN_GEOM = '''    annotate("label", x = x_ann, y = y_ann, label = ann,
             hjust = 0, vjust = 1, size = 2.6,
             fill = "white", alpha = 0.85,
             color = RICE_BLUE_DARK, family = "mono") +'''

# ------------------------------------------------------------------
# AR1 panel (already has xy_lim arg, fix coord block + drop label.size)
# ------------------------------------------------------------------
old_ar1 = '''  # Annotation placed in the upper-left corner using the shared limits
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
             color = RICE_BLUE_DARK, family = "mono") +'''

new_ar1 = NEW_ANNOTATE + '''
  p <- ggplot(d, aes(x = x, y = y)) +
    geom_point(alpha = 0.15, size = 0.5, color = RICE_BLUE) +
    geom_abline(slope = 1, intercept = 0, color = HHD_ORANGE,
                linetype = "dashed", linewidth = 0.5) +
''' + NEW_ANN_GEOM

# ------------------------------------------------------------------
# AR2: full make_scatter_panel rewrite (add xy_lim arg + coord_fixed)
# ------------------------------------------------------------------
# AR2 currently has the rename + annotation but no xy_lim arg.
# Find and replace the entire function body.
ar2_old_body_start = '''make_scatter_panel <- function(df, a, b, log_log = FALSE, title_extra = "") {'''
ar2_new_body_start = '''make_scatter_panel <- function(df, a, b, log_log = FALSE, title_extra = "",
                                xy_lim = NULL) {'''

# Also fix the coord_equal -> coord_fixed switch and annotation block
# AR2's current pattern (from earlier edits): has coord_equal at the end
ar2_old_coord = '''    coord_equal() +
    labs(
      title    = NULL,
      subtitle = NULL,
      x = x_lab, y = y_lab
    ) +'''
ar2_new_coord = '''    labs(
      title    = NULL,
      subtitle = NULL,
      x = x_lab, y = y_lab
    ) +'''

ar2_old_tail = '''  if (log_log) p <- p + scale_x_log10() + scale_y_log10()
  p
}'''
ar2_new_tail = '''  if (log_log) {
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

# AR2 annotation block from earlier rename script
ar2_old_ann = '''  x_ann <- min(d$x, na.rm = TRUE)
  y_ann <- max(d$y, na.rm = TRUE)
  p <- ggplot(d, aes(x = x, y = y)) +
    geom_point(alpha = 0.15, size = 0.5, color = RICE_BLUE) +
    geom_abline(slope = 1, intercept = 0, color = HHD_ORANGE,
                linetype = "dashed", linewidth = 0.5) +
    annotate("label", x = x_ann, y = y_ann, label = ann,
             hjust = 0, vjust = 1, size = 2.6,
             label.size = 0, fill = "white", alpha = 0.85,
             color = RICE_BLUE_DARK, family = "mono") +'''

ar2_new_ann = NEW_ANNOTATE + '''
  p <- ggplot(d, aes(x = x, y = y)) +
    geom_point(alpha = 0.15, size = 0.5, color = RICE_BLUE) +
    geom_abline(slope = 1, intercept = 0, color = HHD_ORANGE,
                linetype = "dashed", linewidth = 0.5) +
''' + NEW_ANN_GEOM

# Also AR2 needs the call sites updated to compute and pass xy_lim
ar2_old_rt_call = '''rt_panels <- lapply(combos, function(cb) {
  make_scatter_panel(rt_mean, cb[1], cb[2], log_log = FALSE,
                     title_extra = "R_t posterior median",
                     xy_lim = rt_xy_lim)
})'''
ar2_new_rt_call = '''rt_all <- unlist(lapply(combos, function(cb) {
  c(rt_mean[[cb[1]]], rt_mean[[cb[2]]])
}))
rt_xy_lim <- c(0, ceiling(max(rt_all, na.rm = TRUE) * 10) / 10)
rt_panels <- lapply(combos, function(cb) {
  make_scatter_panel(rt_mean, cb[1], cb[2], log_log = FALSE,
                     title_extra = "R_t posterior median",
                     xy_lim = rt_xy_lim)
})'''

ar2_old_it_call = '''it_panels <- lapply(combos, function(cb) {
  make_scatter_panel(it_mean, cb[1], cb[2], log_log = TRUE,
                     title_extra = "I_t posterior median, log10",
                     xy_lim = it_xy_lim)
})'''
ar2_new_it_call = '''it_all <- unlist(lapply(combos, function(cb) {
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

# ------------------------------------------------------------------
# Apply patches
# ------------------------------------------------------------------
for f in files:
    p = Path(f)
    s = p.read_text()
    is_ar2 = "_AR2" in f or "ar2" in f.lower()
    orig = s

    # AR1 panel fix (idempotent on AR1)
    s = s.replace(old_ar1, new_ar1)

    # AR2 specific
    if is_ar2:
        # Add xy_lim arg to function signature
        if "xy_lim = NULL" not in s:
            s = s.replace(ar2_old_body_start, ar2_new_body_start)
        # Replace annotation coord block + geom (idempotent if already patched)
        s = s.replace(ar2_old_ann, ar2_new_ann)
        # Remove coord_equal() in the labs block
        s = s.replace(ar2_old_coord, ar2_new_coord)
        # Replace tail to add coord_fixed conditionals
        if "coord_fixed(xlim" not in s:
            s = s.replace(ar2_old_tail, ar2_new_tail)
        # rt/it call sites get the global limits computed
        # (only update if not already computed)
        if "rt_xy_lim <-" not in s:
            s = s.replace(ar2_old_rt_call, ar2_new_rt_call)
        if "it_xy_lim <-" not in s:
            s = s.replace(ar2_old_it_call, ar2_new_it_call)

    p.write_text(s)
    name = f.split('/')[-1]
    if s != orig:
        print(f"  patched {name}")
    else:
        print(f"  NO CHANGE in {name} (patterns may have shifted)")

print()
print("Done. Re-run the two analyses:")
print('  Rscript "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/analyze_models_A_AR1.R"')
print('  Rscript "Modelos AR con rho est. y rho=0/paper_pipeline/08_analysis_AR2/analyze_models_A_AR2.R"')
PYEOF
