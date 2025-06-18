#!/bin/bash
# ----------------SLURM Parameters----------------
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=50GB
#SBATCH --time=15:00:00
#SBATCH --partition=open
#SBATCH -J WGCNA2a2

# ----------------Load Modules--------------------
module load r

# ----------------Commands------------------------
Rscript WGCNAscript2a2.R $1 $2 $3 $4 $5 $6
