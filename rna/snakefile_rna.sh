#!/bin/bash
#SBATCH --job-name=snakemake_rna    # Nombre del trabajo
#SBATCH --output=RNA_%j.out         # Salida
#SBATCH --error=RNA_%j.err          # Errores
#SBATCH --time=10:00:00             # Tiempo máximo de ejecución (HH:MM:SS)
#SBATCH --nodes=1                   # Número de nodos
#SBATCH --ntasks=4                  # Número de tareas (procesos)
#SBATCH --cpus-per-task=4           # CPUs por tarea
#SBATCH --mem=32G                   # Memoria por nodo

module load cesga/2020 snakemake/7.32.4-python-3.10.8
snakemake --slurm --jobs 4 --cores 512 --rerun-incomplete