#!/bin/bash
#SBATCH --job-name=parse_antismash
#SBATCH --output=parse_antismash.out
#SBATCH --error=parse_antismash.err
#SBATCH --time=48:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=10

module load python/3.9.16

INPUT_DIR=/path/to/antismash_output
OUTPUT_CSV=bgc_summary_unique_ids_v8.csv

python 01_parse_antismash_gbk.py \
  --input_dir "${INPUT_DIR}" \
  --output_csv "${OUTPUT_CSV}"
