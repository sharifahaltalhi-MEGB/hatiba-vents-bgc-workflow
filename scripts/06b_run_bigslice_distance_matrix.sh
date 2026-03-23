#!/bin/bash
#SBATCH --job-name=calculate_distance
#SBATCH --output=calculate_distance.%j.out
#SBATCH --error=calculate_distance.%j.err
#SBATCH --time=70:00:00
#SBATCH --mem=900G
#SBATCH --cpus-per-task=8

set -euo pipefail
module purge

# Activate environment
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate bigslice_new

export PYTHONNOUSERSITE=1

# Run distance calculation
python 06_calculate_bigslice_distance_matrix.py
