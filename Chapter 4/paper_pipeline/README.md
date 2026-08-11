# Paper pipeline — Houston RSV spatial Bayesian model

Self-contained pipeline for the paper *"Inferring Spatial Transmission Dynamics
of Respiratory Syncytial Virus Across Houston Wastewater Treatment Plants"*
(Palacio et al., 2026). Two parallel sub-pipelines coexist here:

1. **PEHD pipeline** — the paper's primary spatial weight matrix, built from
   the Population-Extended Hausdorff dissimilarity.
2. **Euclidean pipeline** — an ablation baseline using centroid-to-centroid
   Euclidean distance through the identical RBF kernel + row-normalization
   construction. Used to defend the population-weighting choice against
   reviewer pushback.

Both pipelines share Stan model code, data preprocessing, and HMC settings;
only the spatial weight matrix `W` differs.

## Folder structure

```
paper_pipeline/
├── 01_W_matrix/
│   └── compute_W_matrices.R         # Builds W_PEHD and W_eucl from polygons
├── 02_fitting_PEHD/                 # Primary pipeline (PEHD W)
│   ├── run_A_fsv_est.R              # Joint AR(1) phi-est
│   ├── run_A_marss_est.R            # SSM   AR(1) phi-est
│   ├── run_A_fsv_est.slurm
│   └── run_A_marss_est.slurm
├── 03_fitting_eucl/                 # Ablation pipeline (Euclidean W)
│   ├── run_A_fsv_est_eucl.R         # Joint AR(1) phi-est, Euclidean W
│   ├── run_A_marss_est_eucl.R       # SSM   AR(1) phi-est, Euclidean W
│   ├── run_A_fsv_est_eucl.slurm
│   └── run_A_marss_est_eucl.slurm
├── 04_analysis/
│   ├── analyze_models_A_AR1.R       # CRPS for the 4 PEHD AR(1) fits
│   ├── analyze_models_A_AR1_eucl.R  # CRPS for the 4 Euclidean AR(1) fits
│   ├── analyze_episewer_rw_vs_joint.R   # Joint (PEHD) vs EpiSewer baseline
│   └── analyze_pehd_vs_eucl.R       # Side-by-side PEHD vs Euclidean
├── 05_episewer_baseline/
│   ├── run_episewer_rw.R            # EpiSewer fit, one plant per call
│   ├── run_episewer_rw_treedepth20.R# Fallback for hard plants
│   └── check_convergence.R          # Diagnostics for the 15 EpiSewer fits
└── results/
    ├── W_PEHD.rds                   # Saved by 01_W_matrix
    ├── W_eucl.rds                   # Saved by 01_W_matrix
    ├── d_PEHD.rds, d_eucl.rds       # Distance matrices (km)
    ├── pehd_vs_eucl_crps.csv        # Saved by analyze_pehd_vs_eucl.R
    ├── pehd_vs_eucl_phi.csv         # Saved by analyze_pehd_vs_eucl.R
    └── fig_pehd_vs_eucl_crps.pdf
```

## Workflow

### Step 1 — Build spatial weight matrices (once)

```bash
cd 01_W_matrix
Rscript compute_W_matrices.R
```

Produces `../results/{W_PEHD.rds, W_eucl.rds, d_PEHD.rds, d_eucl.rds}`.

### Step 2a — Fit the PEHD pipeline (primary)

```bash
cd ../02_fitting_PEHD
Rscript run_A_fsv_est.R   2 2 2000 8000   # Joint phi-est
Rscript run_A_marss_est.R 2 2 2000 8000   # SSM   phi-est
# (phi=0 variants come from older v3 fits — see note below)
```

Outputs land in `outputs/stan_fits_fsv_rhoest/` and `outputs/stan_fits_marss_rhoest/`.

### Step 2b — Fit the Euclidean ablation pipeline

```bash
cd ../03_fitting_eucl
Rscript run_A_fsv_est_eucl.R   2 2 2000 8000
Rscript run_A_marss_est_eucl.R 2 2 2000 8000
```

Outputs land in `outputs/stan_fits_fsv_rhoest_eucl/` and `outputs/stan_fits_marss_rhoest_eucl/`.

### Step 3 — Analyze CRPS per pipeline

```bash
cd ../04_analysis
Rscript analyze_models_A_AR1.R         # Writes results_modelA_AR1/
Rscript analyze_models_A_AR1_eucl.R    # Writes results_modelA_AR1_eucl/
```

### Step 4 — EpiSewer external baseline (PEHD-pipeline only)

```bash
cd ../05_episewer_baseline
for IDX in $(seq 1 15); do
  Rscript run_episewer_rw.R $IDX 3000 3000
done
Rscript check_convergence.R
cd ../04_analysis
Rscript analyze_episewer_rw_vs_joint.R
```

### Step 5 — PEHD vs Euclidean side-by-side comparison

```bash
cd 04_analysis
Rscript analyze_pehd_vs_eucl.R
```

Outputs to `../results/`:
- `pehd_vs_eucl_crps.csv` — per-plant comparison table
- `pehd_vs_eucl_phi.csv` — phi posteriors, both pipelines
- `fig_pehd_vs_eucl_crps.pdf` — bar chart

## Notes on the phi=0 (fix) variants

The paper reports four AR(1) variants: `{Joint, SSM} × {phi-est, phi=0}`.
However, the current R driver scripts (`run_A_fsv_est.R`, `run_A_marss_est.R`)
only fit the phi-estimated versions. The phi=0 fits in the paper are
**inherited from an earlier pipeline run (the "v3" series)**:

- `fits/A_joint_fix.rds` → symlink to `outputs/stan_fits_fsv_v3/rho_00.rds`
- `fits/A_ssm_fix.rds`   → symlink to `outputs/stan_fits_ensor_v3/rho_00.rds`

For the Euclidean ablation to include all four variants, two additional driver
scripts (`run_A_fsv_fix_eucl.R`, `run_A_marss_fix_eucl.R`) would be needed.
These would use the Stan models `renewal_A_fsv.stan` and `renewal_A_marss.stan`
(which take `rho` as data rather than a parameter) with `rho = 0`. Decide
whether the Euclidean ablation needs the phi=0 variants or only the
phi-estimated ones, since the paper's primary CRPS comparison
(Joint phi-est vs EpiSewer) only uses phi-est.

## Reproducing the headline CRPS result

After fitting all four PEHD variants and the 15 EpiSewer baselines:

```bash
cd 04_analysis
Rscript analyze_models_A_AR1.R         # produces crps_by_plant_AR1.csv
Rscript analyze_episewer_rw_vs_joint.R # uses it for the final comparison
```

Final aggregate: **Joint AR(1) with phi-est outperforms EpiSewer in 15 of 15
plants, mean log10 CRPS reduction = 30.6% (range 22.1–46.3%).**

## R package requirements

```r
install.packages(c(
  "tidyverse", "lubridate", "sf", "MARSS",
  "cmdstanr", "posterior",
  "scoringRules", "loo",
  "ggplot2", "ggh4x", "scales", "gridExtra"
))
remotes::install_github("adrian-lison/EpiSewer")
```

`cmdstanr` requires a working CmdStan installation.

## Contact

Jose R. Palacio (`jrp16@rice.edu`) — Department of Statistics, Rice University.
