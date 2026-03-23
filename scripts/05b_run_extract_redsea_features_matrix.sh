#!/bin/bash
#SBATCH --job-name=extract_bigslice_features
#SBATCH --output=extract_bigslice_features.%j.out
#SBATCH --error=extract_bigslice_features.%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=2

set -euo pipefail
module purge

# Activate environment
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate bigslice_new

export PYTHONNOUSERSITE=1

# Input: BiG-SLiCE result directory (e.g. run folder)
BIGSLICE_RESULT_DIR="/path/to/bigslice_run_directory"

# Output: feature matrix
OUTPUT_MATRIX="RedSea_features_matrix.tsv"

# Run feature extraction
python 05_extract_bigslice_features_matrix.py \
  "$BIGSLICE_RESULT_DIR" \
  "$OUTPUT_MATRIX"
