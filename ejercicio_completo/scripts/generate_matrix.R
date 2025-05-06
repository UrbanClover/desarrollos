.libPaths(c("/mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/R/library", .libPaths()))

# Cargar librerías necesarias
library(tximport)
library(readr)
library(dplyr)
library(tibble)
library(stringr)

# Accede a inputs y outputs desde Snakemake
input_files <- snakemake@input[1:(length(snakemake@input)-1)]
metadata_file <- snakemake@input[[length(snakemake@input)]]
output_file <- snakemake@output[[1]]

# Configuración de directorios
salmon_dir <- "results/salmon"
output_dir <- "results/expression_matrix"
samples <- list.dirs(salmon_dir, recursive = FALSE, full.names = FALSE)

# Obtener rutas de los archivos quant.sf
files <- file.path(salmon_dir, samples, "quant.sf")
names(files) <- samples

# Importar datos de Salmon
txi <- tximport(files, type = "salmon", txOut = TRUE)

# Crear matriz de expresión (TPM)
expression_matrix <- as.data.frame(txi$counts) %>%
  rownames_to_column("transcript") %>%
  mutate(across(where(is.numeric), round, 2))

# Guardar matriz
dir.create(output_dir, showWarnings = FALSE)
write_tsv(expression_matrix, file.path(output_dir, "transcript_matrix.tsv"))

message("Matriz de expresión generada exitosamente!")