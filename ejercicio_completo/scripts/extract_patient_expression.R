#!/usr/bin/env Rscript

# ==============================================================================
# Script: extract_patient_expression.R (Versión sin errores de sintaxis)
# ==============================================================================

# --- 1. Configuración inicial ---
library(optparse)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

# --- 2. Manejo de argumentos ---
option_list <- list(
  make_option(c("-m", "--matrix"), type = "character", help = "Matriz de expresión (transcript_matrix.tsv)"),
  make_option(c("-t", "--top"), type = "character", help = "Top transcritos (top10_common_transcripts.txt)"),
  make_option(c("-s", "--metadata"), type = "character", help = "Metadatos (sample_metadata.tsv)"),
  make_option(c("-o", "--output"), type = "character", help = "Archivo de salida (patient_expression_values.tsv)")
)
opt <- parse_args(OptionParser(option_list = option_list))

# --- 3. Carga de datos ---
expression_data <- fread(opt$matrix) %>%
  pivot_longer(
    cols = -transcript,
    names_to = "sample",
    values_to = "tpm"
  ) %>%
  mutate(
    transcript = str_split(transcript, "\\|", simplify = TRUE)[,1]  # Normalizar nombres
  )

# Leer y procesar top_transcripts (corregido)
top_transcripts <- readLines(opt$top) %>%
  str_split("\\|") %>%
  sapply(function(x) x[1])  # Función anónima explícita

metadata <- fread(opt$metadata)

# --- 4. Procesamiento de datos ---
processed_data <- expression_data %>%
  mutate(
    patient_id = str_extract(sample, "^\\d+"),  # Extraer "001", "002", etc.
    tumor_type = str_remove(sample, "^\\d+_")   # Ej: "primary_tumor"
  ) %>%
  filter(transcript %in% top_transcripts) %>%
  group_by(patient_id, tumor_type, transcript) %>%
  summarise(
    mean_tpm = round(mean(tpm, na.rm = TRUE)), 
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = transcript,
    values_from = mean_tpm,
    values_fill = list(mean_tpm = 0)
  )

# --- 5. Integración con metadatos ---
final_output <- metadata %>%
  separate(
    sample_id, 
    into = c("patient_id", "tumor_type"), 
    sep = "_", 
    extra = "merge"
  ) %>%
  left_join(processed_data, by = c("patient_id", "tumor_type"))

# --- 6. Guardar resultados ---
fwrite(final_output, opt$output, sep = "\t")