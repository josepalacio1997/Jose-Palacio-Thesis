#!/bin/bash
#SBATCH --job-name=insample_rhoest_ar2
#SBATCH --partition=commons
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-user=jrp16@nots.rice.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --no-requeue

set -euo pipefail
echo "=== SLURM job $SLURM_JOB_ID ($SLURM_JOB_NAME) start: $(date) ==="
echo "Host: $(hostname)"
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs outputs/stan_fits_marss_rhoest_ar2

module purge
module load GCC/13.2.0
module load OpenMPI/4.1.6
module load R/4.4.1
module load R-bundle-CRAN/2024.06 2>/dev/null || true
module load R-bundle-Bioconductor/3.19-R-4.4.1 2>/dev/null || true
module load GDAL/3.9.0 GEOS/3.12.1 PROJ/9.3.1

echo "R version: $(Rscript --version 2>&1)"
echo "Free mem:"; free -h

# MARSS AR(2) con rho ESTIMADO, prior Beta(2,2), 2000 warmup + 8000 sampling
# Incluye MARSS Stage 1 + Stan rho-est AR(2) con (phi1, phi2) = (2, -1)
Rscript --no-save --no-restore run_A_marss_est_ar2.R 2 2 2000 8000

echo "=== SLURM job $SLURM_JOB_ID end: $(date) ==="
ls -la outputs/stan_fits_marss_rhoest_ar2/rho_estimated.rds 2>/dev/null || true
