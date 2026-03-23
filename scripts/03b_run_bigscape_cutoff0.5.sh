#!/bin/bash
#SBATCH --job-name=bigscape_c0p7
#SBATCH --output=bigscape_c0p7.out
#SBATCH --error=bigscape_c0p7.err
#SBATCH --time=150:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=40

set -eo pipefail

source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate bigscape2_env

INPUT_DIR="/path/to/gbk_by_genome"
OUTPUT_DIR="/path/to/bigscape_cutoff0p7_output"
PFAM_PATH="/path/to/Pfam-A.hmm"

[[ -d "$INPUT_DIR" ]] || { echo "ERROR: INPUT_DIR not found: $INPUT_DIR"; exit 1; }
[[ -f "$PFAM_PATH" ]] || { echo "ERROR: PFAM_PATH not found: $PFAM_PATH"; exit 1; }

bigscape --version

bigscape cluster \
  -i "$INPUT_DIR" \
  -o "$OUTPUT_DIR" \
  -c 40 \
  --input-mode recursive \
  --mix \
  --include-singletons \
  --gcf-cutoffs 0.5 \
  --alignment-mode glocal \
  -p "$PFAM_PATH" \
  -m 3.1 \
  -v

echo "BiG-SCAPE 2 complete: $OUTPUT_DIR"
