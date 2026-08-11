#!/usr/bin/env bash
# =============================================================================
# drop_rmse_from_scatter.sh
# -----------------------------------------------------------------------------
# Remove RMSE from the scatter-panel annotations, leaving only the Pearson r.
# Matches the simplified figure captions in the paper.
#
# Run from the repo root:
#   bash "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/drop_rmse_from_scatter.sh"
# =============================================================================

set -e
cd "$HOME/Documents/Claude/Projects/NOTS/rsv_mcmc_stan_geoW"

FILES=(
  "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/analyze_models_A_AR1.R"
  "Modelos AR con rho est. y rho=0/paper_pipeline/08_analysis_AR2/analyze_models_A_AR2.R"
)

DATE=$(date +%Y%m%d_droprmse)
for f in "${FILES[@]}"; do
  cp "$f" "$f.bak.$DATE"
done

python3 << 'PYEOF'
from pathlib import Path

files = [
    "Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/analyze_models_A_AR1.R",
    "Modelos AR con rho est. y rho=0/paper_pipeline/08_analysis_AR2/analyze_models_A_AR2.R",
]

# Replace the multi-line annotation block with a single-line "r" only
old1 = '''  r    <- cor(xx, yy)
  rmse <- sqrt(mean((xx - yy)^2))
  ann  <- sprintf("r = %.3f\\nRMSE = %.3f", r, rmse)'''
new1 = '''  r   <- cor(xx, yy)
  ann <- sprintf("r = %.3f", r)'''

for f in files:
    p = Path(f)
    s = p.read_text()
    if old1 in s:
        s = s.replace(old1, new1)
        p.write_text(s)
        print(f"  patched {f.split('/')[-1]}")
    else:
        print(f"  NOT MATCHED in {f.split('/')[-1]}")

print("\nDone. Re-run the analyses to regenerate figures with only r:")
print("  Rscript \"Modelos AR con rho est. y rho=0/paper_pipeline/04_analysis/analyze_models_A_AR1.R\"")
print("  Rscript \"Modelos AR con rho est. y rho=0/paper_pipeline/08_analysis_AR2/analyze_models_A_AR2.R\"")
PYEOF
