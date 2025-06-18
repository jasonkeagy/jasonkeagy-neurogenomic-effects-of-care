#!/bin/bash
# ----------------SLURM Parameters----------------
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=5GB
#SBATCH --time=1:00:00
#SBATCH --partition=open
#SBATCH -J WGCNA2a1

# ----------------Load Modules--------------------
module load r

# ----------------Commands------------------------

Rscript WGCNAscript2a1.R $1
