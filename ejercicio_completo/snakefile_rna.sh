#!/bin/bash
#SBATCH --job-name=Snakefile_rna
#SBATCH --output=Snakefile_rna_%j.out
#SBATCH --error=Snakefile_rna_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --account=othcx

module load cesga/2020 snakemake/7.32.4-python-3.10.8
snakemake -s Snakefile_rna "$@" --slurm --jobs 10 --cores 254 --rerun-incomplete