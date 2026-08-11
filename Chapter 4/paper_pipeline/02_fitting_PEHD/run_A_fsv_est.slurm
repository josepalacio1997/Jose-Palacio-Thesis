#!/bin/bash
#SBATCH --job-name=insample_fsv_rhoest
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
mkdir -p logs outputs/stan_fits_fsv_rhoest

module purge
module load GCC/13.2.0
module load OpenMPI/4.1.6
module load R/4.4.1
module load R-bundle-CRAN/2024.06 2>/dev/null || true
module load R-bundle-Bioconductor/3.19-R-4.4.1 2>/dev/null || true
module load GDAL/3.9.0 GEOS/3.12.1 PROJ/9.3.1

echo "R version: $(Rscript --version 2>&1)"
echo "Free mem:"; free -h

# FSV con rho ESTIMADO, prior Beta(2,2), 2000 warmup + 8000 sampling
# (mismo config MCMC que los fits rho fijo, para comparabilidad)
Rscript --no-save --no-restore run_A_fsv_est.R 2 2 2000 8000

echo "=== SLURM job $SLURM_JOB_ID end: $(date) ==="
ls -la outputs/stan_fits_fsv_rhoest/rho_estimated.rds 2>/dev/null || true
