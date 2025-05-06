# --------------------------------------
# 1. Configuración inicial
# --------------------------------------
.libPaths(c("/mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/R/library", .libPaths()))
# Cargar paquetes
library(tidyverse)
library(ggpubr)
library(gridExtra)

# Definir muestras y rutas
muestras <- paste0("LC_", 1:5)  # LC_1, LC_2, ..., LC_5
directorio_base <- "resultados/annovar/strelka"

# Obtener rutas de los archivos multianno.txt
files <- file.path(
  directorio_base,
  muestras,
  paste0(muestras, ".hg38_multianno.txt")
)

# Verificar que todos los archivos existen
if(!all(file.exists(files))) {
  stop("¡Algunos archivos no se encuentran! Verifica las rutas.")
}

# --------------------------------------
# 2. Leer y combinar todos los archivos
# --------------------------------------
# Función para leer cada archivo con información de la muestra
leer_annovar <- function(file) {
  muestra <- strsplit(file, "/")[[1]][4]  # Extraer LC_X del path
  df <- read_delim(file, delim = "\t", show_col_types = FALSE) %>%
    mutate(Muestra = muestra)  # Añadir columna de identificación
  return(df)
}

# Leer todos los archivos y combinar en un solo dataframe
df_completo <- map_dfr(files, leer_annovar)

# Verificar que las columnas necesarias existen
required_columns <- c("ExonicFunc.refGene", "Func.refGene", "Chr", "Start", "Ref", "Alt")
missing_columns <- setdiff(required_columns, colnames(df_completo))
if (length(missing_columns) > 0) {
  stop(paste("Las siguientes columnas faltan en el dataframe:", paste(missing_columns, collapse = ", ")))
}

# --------------------------------------
# 3. Análisis conjunto
# --------------------------------------
# A. Resumen general
resumen_general <- df_completo %>%
  group_by(Muestra) %>%
  summarise(
    Total_Variantes = n(),
    Variantes_Exonicas = sum(Func.refGene == "exonic", na.rm = TRUE),
    Porcentaje_Exonicas = round(100 * Variantes_Exonicas / Total_Variantes, 1)
  )

print(resumen_general)

# B. Gráfico combinado de distribución de variantes
grafico_conjunto <- ggplot(df_completo, aes(x = Func.refGene, fill = Muestra)) +
  geom_bar(position = "dodge") +
  labs(
    title = "Distribución de Variantes por Región Génica",
    x = "Región Génica",
    y = "Número de Variantes"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(grafico_conjunto)

# C. Variantes compartidas entre muestras
variantes_unicas <- df_completo %>%
  unite("ID_Variante", Chr, Start, Ref, Alt, sep = "-") %>%
  group_by(ID_Variante) %>%
  summarise(
    Muestras = paste(unique(Muestra), collapse = ", "),
    Conteo = n()
  ) %>%
  arrange(desc(Conteo))

# Top 10 variantes más frecuentes
print(head(variantes_unicas, 10))

# --------------------------------------
# 4. Análisis individual por muestra
# --------------------------------------
# Crear lista de gráficos para cada muestra
graficos_individuales <- df_completo %>%
  group_by(Muestra) %>%
  group_map(~ {
    ggplot(.x, aes(x = ExonicFunc.refGene, fill = ExonicFunc.refGene)) +
      geom_bar() +
      labs(
        title = paste("Variantes Exónicas -", .y$Muestra),
        x = "Tipo de Variante",
        y = "Conteo"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })

# Combinar todos los gráficos en un PDF
ggsave(
  "Analisis_Individual_Muestras.pdf",
  gridExtra::marrangeGrob(graficos_individuales, nrow = 2, ncol = 2),  
  width = 15,
  height = 9
)

# --------------------------------------
# 5. Exportar datos combinados
# --------------------------------------
write_csv(df_completo, "Variantes_Combinadas_LC.csv")