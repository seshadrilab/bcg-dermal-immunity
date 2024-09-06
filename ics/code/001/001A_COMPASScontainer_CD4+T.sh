#!/bin/bash
#SBATCH --partition=campus-new
#SBATCH --cpus-per-task=32
#SBATCH --time=72:00:00
#SBATCH --mem=300000
#SBATCH -o ./log_001A_COMPASScontainer_CD4+T
#SBATCH -J BCG-TICE

echo "Start of program at `date`"

ml fhR/4.1.0-foss-2020b
Rscript --no-save --no-restore ./001A_COMPASScontainer_CD4+T.R

echo "End of program at `date`"
