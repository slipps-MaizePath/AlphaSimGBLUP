#!/usr/bin/env bash
#
#SBATCH --job-name=SimulationsBig
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=65:00:00
#SBATCH --partition=IllinoisComputes
#SBATCH --account=slipps-ic
#SBATCH --array=1-593%100


module load R/4.5.1-tidyverse

Rscript ./02_RunPipeline.R
