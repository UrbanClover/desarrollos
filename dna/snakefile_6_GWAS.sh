#!/bin/bash
#SBATCH --job-name=Snakefile_6_GWAS
#SBATCH --output=Snakefile_6_GWAS_%j.out
#SBATCH --error=Snakefile_6_GWAS_%j.err
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --account=othcx

module purge
module load cesga/2020 snakemake/7.32.4-python-3.10.8 gcc/system plink/2.00a2.3

snakemake -s Snakefile_6_GWAS -p --jobs 5 --cores 4 --latency-wait 30