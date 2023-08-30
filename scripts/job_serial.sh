#!/bin/bash
#SBATCH --job-name=donors
#SBATCH --output=donors_%j.out
#SBATCH --error=donors_%j.err

#SBATCH --ntasks=1
#SBATCH --qos=bsc_ls
#SBATCH --time=2-00:00:00

module load gcc pcre2 R/4.3.0
module load python/3.6.1

mkdir ../outputs
Rscript 00_parse_tree.R
python3 01_get_branch_set.py
Rscript 02_simulate_transfers.R
Rscript 03_plotting_shifts.R
Rscript 04_groups_data.R
Rscript 05_heatmap.R

