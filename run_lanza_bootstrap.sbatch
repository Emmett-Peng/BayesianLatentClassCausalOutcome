#!/bin/bash
#SBATCH --job-name=lanza_boot
#SBATCH --account=def-YOUR_ACCOUNT   # replace with your allocation account
#SBATCH --array=1-1600%200           # 16 scenarios x 100 datasets; 200 tasks at a time
#SBATCH --cpus-per-task=1            # one dataset per task; flexmix is single-core
#SBATCH --mem=4G
#SBATCH --time=24:00:00              # about 4 hours per dataset at B = 1000
#SBATCH --output=logs/lanza_%A_%a.out
#SBATCH --error=logs/lanza_%A_%a.err

# Comparator (two-step IPW latent class estimator) with B = 1000 bootstrap replicates.
# Each array task handles one (scenario, dataset) pair and writes its own checkpoint,
# so a failed or timed-out task can simply be resubmitted: finished datasets are skipped.
#
#   mkdir -p logs && sbatch run_lanza_bootstrap.sbatch
#
# If the cluster refuses a 1600-element array, submit it in chunks (the task-to-dataset
# mapping is fixed, so chunks are safe):
#   sbatch --array=1-800%200    run_lanza_bootstrap.sbatch
#   sbatch --array=801-1600%200 run_lanza_bootstrap.sbatch
#
# When all tasks have finished, build the two CSV files:
#   module load r/4.5.0 && AGGREGATE_ONLY=1 Rscript lanza_bootstrap_B1000.R
#
# B is set in lanza_bootstrap_B1000.R (B_BOOT, default 1000), so there is one place to change it.

module load r/4.5.0                  # adjust to the R module available on your cluster

set -eo pipefail
if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript not on PATH. Check the 'module load' line above." >&2
  exit 1
fi

export WORKERS=1

Rscript lanza_bootstrap_B1000.R
