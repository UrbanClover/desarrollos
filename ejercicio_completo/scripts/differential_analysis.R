#!/usr/bin/env Rscript

# ====================================================================================
# Script: differential_analysis.R
# Descripción: Análisis de expresión diferencial usando DESeq2 a partir de matrices
#              generadas por Snakemake. 
# Uso: Rscript differential_analysis.R -i transcript_matrix.tsv -m metadata.tsv -o top10.txt -d "~tumor_type"
# ====================================================================================

# --- 1. Configuración inicial ---
library(optparse)
library(tidyverse)
library(DESeq2)

# Establecer directorio de librerías personalizado (¡No instalar paquetes durante la ejecución!)
.libPaths(c("/mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/R/library", .libPaths()))

# --- 2. Manejo de argumentos ---
option_list <- list(
  make_option(c("-i","--input"),  type = "character", help = "Matriz de expresión (transcript_matrix.tsv)"),
  make_option(c("-m","--meta"),   type = "character", help = "Metadatos (sample_metadata.tsv)"),
  make_option(c("-o","--output"), type = "character", help = "Archivo de salida con top10 transcritos"),
  make_option(c("-d","--design"), type = "character", help = "Fórmula de diseño (e.g., '~ tumor_type')")
)
opt <- parse_args(OptionParser(option_list = option_list))

# --- 3. Carga de datos ---
# Matriz de expresión (transcriptos x muestras)
count_matrix <- read_tsv(opt$input) %>%
  column_to_rownames("transcript") %>%  # Los transcritos deben ser los nombres de fila
  as.matrix()

# Metadatos (muestras x variables)
metadata <- read_tsv(opt$meta) %>%
  mutate(tumor_type = factor(tumor_type))  # Convertir a factor la variable de interés

# --- 4. Verificación de integridad ---
# Asegurar que los nombres de las columnas de la matriz coincidan con los metadatos
stopifnot(all(colnames(count_matrix) %in% metadata$sample_id))

# --- 5. Análisis con DESeq2 ---
# Crear objeto DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),  # DESeq2 requiere valores enteros
  colData   = metadata,
  design    = as.formula(opt$design)
)

# Filtrado de transcritos de baja expresión (¡ajusta según tus datos!)
keep <- rowSums(counts(dds) >= 10) >= 3  # Conservar transcritos con ≥10 cuentas en ≥3 muestras
dds <- dds[keep,]

# Normalización y análisis diferencial
dds <- DESeq(dds)

# Extraer resultados (primary_tumor vs recurrent_tumor)
res <- results(
  dds,
  contrast = c("tumor_type", "primary_tumor", "recurrent_tumor"),
  alpha = 0.05
)

# --- 6. Selección de transcritos significativos ---
top_transcripts <- res %>%
  as.data.frame() %>%
  rownames_to_column("transcript") %>%
  arrange(padj) %>%  # Ordenar por p-valor ajustado
  filter(!is.na(padj)) %>%  # Eliminar NA
  head(10) %>%  # Top 10 transcritos
  pull(transcript)

# --- 7. Guardar resultados ---
write_lines(top_transcripts, opt$output)