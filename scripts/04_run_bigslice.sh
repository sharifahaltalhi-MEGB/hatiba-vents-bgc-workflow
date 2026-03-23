#!/bin/bash
#SBATCH --job-name=bigslice
#SBATCH --output=bigslice.%j.out
#SBATCH --error=bigslice.%j.err
#SBATCH --time=150:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=40

set -euo pipefail

module purge
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate bigslice_new

export PYTHONNOUSERSITE=1

INPUT_DIR="/path/to/maingbk_by_genome"
BASE_OUT="/path/to/bigslice_output"
OUT_DIR="$BASE_OUT/run_${SLURM_JOB_ID}"

mkdir -p "$BASE_OUT"
rm -rf "$OUT_DIR"

which bigslice
bigslice --version

bigslice -i "$INPUT_DIR" -t "$SLURM_CPUS_PER_TASK" "$OUT_DIR"
