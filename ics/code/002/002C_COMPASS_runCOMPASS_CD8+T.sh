#!/bin/bash
#SBATCH --partition=campus-new
#SBATCH --cpus-per-task=32
#SBATCH --time=148:00:00
#SBATCH --mem=300000
#SBATCH -o ./log_002C_COMPASS_runCOMPASS_CD8+T
#SBATCH -J BCG-TICE

echo "Start of program at `date`"

ml fhR/4.1.0-foss-2020b
Rscript --no-save --no-restore ./002C_COMPASS_runCOMPASS_CD8+T.R

echo "End of program at `date`"
