#!/usr/bin/env bash
# =============================================================================
# rename_variants_and_add_scatter_stats.sh
# -----------------------------------------------------------------------------
# Two changes to the 4 analyze scripts:
#   1. Rename plot labels:  "SSM"   -> "Filter Est"
#                           "Direct" -> "Data"
#      (matches the variant rename Jose applied to main_v4.tex)
#
#   2. Add Pearson r + RMSE annotation to the scatter plot panels
#      (fig9_Rt_scatter_*.pdf and fig10_It_scatter_*.pdf) so the
#      caption's claim that "Pearson correlation and RMSE are reported
#      in each panel" is actually visible in the figures.
#
# Backups (.bak.YYYYMMDD) are created automatically.
#
# Run from the repo root:
#   bash "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/rename_variants_and_add_scatter_stats.sh"
# =============================================================================

set -e

PROJ_ROOT="$HOME/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW"
cd "$PROJ_ROOT"

FILES=(
  "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/analyze_models_A_AR1.R"
  "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/analyze_models_A_AR1_eucl.R"
  "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/analyze_models_A_AR1_HD.R"
  "Modelos AR con rho est. y rho=0/paper_pipeline/08_analysis_AR2/analyze_models_A_AR2.R"
)

# ---- Backup ----
DATE=$(date +%Y%m%d)
for f in "${FILES[@]}"; do
  cp "$f" "$f.bak.$DATE"
done
echo "Backups created with suffix .bak.$DATE"
echo ""

# ---- Apply the patches via Python ----
python3 << 'PYEOF'
import re
from pathlib import Path

ar1_files = [
    "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/analyze_models_A_AR1.R",
    "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/analyze_models_A_AR1_eucl.R",
    "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/analyze_models_A_AR1_HD.R",
]
ar2_file = "Modelos AR con rho est. y rho=0/paper_pipeline/08_analysis_AR2/analyze_models_A_AR2.R"

# ---- Label renames (applied to all 4 scripts) ----
label_replacements = [
    ('bquote("SSM; " * hat(rho))',    'bquote("Filter Est; " * hat(rho))'),
    ('bquote("Direct; " * hat(rho))', 'bquote("Data; " * hat(rho))'),
    ('bquote("SSM; " * rho == 0)',    'bquote("Filter Est; " * rho == 0)'),
    ('bquote("Direct; " * rho == 0)', 'bquote("Data; " * rho == 0)'),
    ('"SSM; rho-est"',    '"Filter Est; rho-est"'),
    ('"Direct; rho-est"', '"Data; rho-est"'),
    ('"SSM; rho=0"',      '"Filter Est; rho=0"'),
    ('"Direct; rho=0"',   '"Data; rho=0"'),
    ('"SSM*\'; \'*hat(rho)"',    '"\'Filter Est; \'*hat(rho)"'),
    ('"Direct*\'; \'*hat(rho)"', '"\'Data; \'*hat(rho)"'),
    ('"SSM*\'; \'*rho==0"',      '"\'Filter Est; \'*rho==0"'),
    ('"Direct*\'; \'*rho==0"',   '"\'Data; \'*rho==0"'),
]

# ---- AR1 scatter panel rewrite ----
old_panel_ar1 = '''make_scatter_panel <- function(df, a, b, log_log = FALSE, title_extra = "") {
  d <- data.frame(x = df[[a]], y = df[[b]]) |>
    filter(is.finite(x) & is.finite(y))
  if (log_log) d <- filter(d, x > 0, y > 0)
  r <- if (log_log) cor(log10(d$x), log10(d$y)) else cor(d$x, d$y)
  # Map fit IDs to formatted axis labels: "SSM; rho-hat" etc.
  x_lab <- title_label[[a]]; if (is.null(x_lab)) x_lab <- a
  y_lab <- title_label[[b]]; if (is.null(y_lab)) y_lab <- b
  p <- ggplot(d, aes(x = x, y = y)) +
    geom_point(alpha = 0.15, size = 0.5, color = RICE_BLUE) +
    geom_abline(slope = 1, intercept = 0, color = HHD_ORANGE,
                linetype = "dashed", linewidth = 0.5) +
    coord_equal() +
    labs(x = x_lab, y = y_lab) +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(face = "bold", color = RICE_BLUE_DARK),
          plot.subtitle = element_text(color = NEUTRAL_GRAY, size = 8))
  if (log_log) p <- p + scale_x_log10() + scale_y_log10()
  p
}'''

new_panel_ar1 = '''make_scatter_panel <- function(df, a, b, log_log = FALSE, title_extra = "") {
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

# ---- Apply to AR1 siblings ----
for f in ar1_files:
    p = Path(f)
    s = p.read_text()
    for old, new in label_replacements:
        s = s.replace(old, new)
    s = s.replace(old_panel_ar1, new_panel_ar1)
    p.write_text(s)
    print(f"  patched {f.split('/')[-1]}")

# ---- AR2: apply label renames (same) ----
p = Path(ar2_file)
s = p.read_text()
for old, new in label_replacements:
    s = s.replace(old, new)
# AR2 has a slightly different make_scatter_panel signature; locate it
# and inject the annotation. Look for the ggplot block.
ar2_old_block = '''  p <- ggplot(d, aes(x = x, y = y)) +
    geom_point(alpha = 0.15, size = 0.5, color = RICE_BLUE) +
    geom_abline(slope = 1, intercept = 0, color = HHD_ORANGE,
                linetype = "dashed", linewidth = 0.5) +
    coord_equal() +
    labs(
      title    = NULL,
      subtitle = NULL,
      x = x_lab, y = y_lab
    ) +'''
ar2_new_block = '''  if (log_log) {
    xx <- log10(d$x); yy <- log10(d$y)
  } else {
    xx <- d$x; yy <- d$y
  }
  r    <- cor(xx, yy)
  rmse <- sqrt(mean((xx - yy)^2))
  ann  <- sprintf("r = %.3f\\nRMSE = %.3f", r, rmse)
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
    labs(
      title    = NULL,
      subtitle = NULL,
      x = x_lab, y = y_lab
    ) +'''
if ar2_old_block in s:
    s = s.replace(ar2_old_block, ar2_new_block)
    print(f"  patched {ar2_file.split('/')[-1]} (with scatter annotation)")
else:
    print(f"  WARN: AR2 scatter block pattern not matched in {ar2_file.split('/')[-1]}")
    print(f"        Label rename applied; scatter annotation may need manual edit.")
p.write_text(s)

print()
print("=== Verification ===")
for f in ar1_files + [ar2_file]:
    txt = Path(f).read_text()
    n_ssm    = txt.count('"SSM;') + txt.count("\"SSM*'")
    n_direct = txt.count('"Direct;') + txt.count("\"Direct*'")
    n_filter = txt.count('Filter Est')
    n_data   = txt.count('"Data;') + txt.count("'Data; '")
    n_rmse   = txt.count("RMSE")
    name = f.split('/')[-1]
    print(f"  {name:40s}  SSM={n_ssm} Direct={n_direct} "
          f"Filter={n_filter} Data={n_data} RMSE-annots={n_rmse}")
PYEOF

echo ""
echo "Done. Now regenerate the figures by re-running the analyses:"
echo ""
echo "  Rscript \"Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/analyze_models_A_AR1.R\""
echo "  Rscript \"Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/analyze_episewer_rw_vs_joint.R\""
echo "  Rscript \"Modelos AR con rho est. y rho=0/paper_pipeline/08_analysis_AR2/analyze_models_A_AR2.R\""
