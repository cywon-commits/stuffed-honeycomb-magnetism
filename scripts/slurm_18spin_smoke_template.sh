#!/bin/bash
#SBATCH --job-name=shc-18smoke
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=shc_18smoke_%j.out
set -euo pipefail
# Adapt Python module/venv and partition/account for your site before sbatch.
cd "$(dirname "$0")/../workspaces/finite_size_18_27_36pilot"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
find . -name '*.so' -delete
python run_all.py --tests
