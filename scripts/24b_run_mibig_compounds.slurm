#!/bin/bash
#SBATCH --job-name=mibig_products
#SBATCH --output=mibig_products_%j.out
#SBATCH --error=mibig_products_%j.err
#SBATCH --time=100:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=1

module load python || true

python 24_get_mibig_compounds.py
