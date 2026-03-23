#!/bin/bash
#SBATCH --job-name=summarize_novelty
#SBATCH --output=summarize_novelty.%j.out
#SBATCH --error=summarize_novelty.%j.err
#SBATCH --time=150:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=2

set -euo pipefail
module purge

source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate bigslice_new
export PYTHONNOUSERSITE=1

python 07_summarize_bigslice_novelty.py
