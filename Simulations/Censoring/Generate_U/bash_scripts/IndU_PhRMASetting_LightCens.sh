#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=5g
#SBATCH -n 3
#SBATCH -t 1-

Rscript IndU_PhRMASetting_LightCens.R

