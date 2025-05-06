#!/bin/bash
#SBATCH --job-name=snakemake_wes    # Nombre del trabajo
#SBATCH --output=snakefile_muestras_normal_%j.out         # Salida
#SBATCH --error=snakefile_muestras_normal_%j.err          # Errores
#SBATCH --time=48:00:00             # Tiempo máximo de ejecución (HH:MM:SS)
#SBATCH --nodes=1                   # Número de nodos
#SBATCH --ntasks=4                  # Número de tareas (procesos)
#SBATCH --cpus-per-task=4           # CPUs por tarea
#SBATCH --mem=4G                   # Memoria por nodo

module load cesga/2020 snakemake/7.32.4-python-3.10.8
snakemake  -s snakefile_muestras_normal "$@" --slurm --jobs 5 --cores 256 --rerun-incomplete