.libPaths(c("/mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/R/R-4.3.0/library", .libPaths()))

library(qqman)
args <- commandArgs(trailingOnly = TRUE)

# ----------------------------
# 1. Leer y filtrar datos
# ----------------------------
gwas_data <- read.table(args[1], header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE)
gwas_data <- subset(gwas_data, TEST == "ADD")
gwas_data$P <- as.numeric(gwas_data$P)
gwas_data <- na.omit(gwas_data)

# Verificar si hay datos
if (nrow(gwas_data) == 0) stop("No hay datos después de filtrar por TEST == ADD y eliminar NA en P.")

# ----------------------------
# 2. Corregir columna X.CHROM
# ----------------------------
# Extraer valores únicos de X.CHROM
unique_chr <- unique(gwas_data$X.CHROM)
print("Valores únicos en X.CHROM antes de conversión:")
print(unique_chr)

# Mapear TODOS los cromosomas no numéricos a números secuenciales
numeric_chr <- suppressWarnings(as.numeric(unique_chr))
non_numeric <- unique_chr[is.na(numeric_chr)]  # Cromosomas no convertibles

# Asignar números únicos a partir de 23
chr_mapping <- setNames(23:(22 + length(non_numeric)), non_numeric)
gwas_data$X.CHROM <- ifelse(
  gwas_data$X.CHROM %in% non_numeric,
  chr_mapping[gwas_data$X.CHROM],
  as.numeric(gwas_data$X.CHROM)  # Cromosomas numéricos originales
)

# Convertir a numérico definitivo y eliminar NAs
gwas_data$X.CHROM <- as.numeric(gwas_data$X.CHROM)
gwas_data <- na.omit(gwas_data)

# Verificar datos post-conversión
print("Valores únicos en X.CHROM después de conversión:")
print(unique(gwas_data$X.CHROM))
if (nrow(gwas_data) == 0) stop("No hay datos válidos después de convertir X.CHROM.")

# ----------------------------
# 3. Generar Manhattan Plot
# ----------------------------
colors <- c("#1E90FF", "#FF4500")
sig_line <- 0.7  # -log10(p) = 0.7 → p = 0.1995

# Asegurar directorio de salida
output_dir <- dirname(args[2])
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

png(args[2], width = 1200, height = 600, res = 150)
manhattan(gwas_data,
          chr = "X.CHROM",
          bp = "POS",
          p = "P",
          snp = "ID",
          col = colors,
          genomewideline = FALSE,
          suggestiveline = sig_line,
          main = "Manhattan Plot",
          cex = 0.6,
          ylim = c(0, max(-log10(gwas_data$P), na.rm = TRUE) * 1.1))

abline(h = sig_line, col = "red", lty = 2)
legend("topleft", legend = paste("Umbral: p <", round(10^-sig_line, 3)), col = "red", lty = 2, bty = "n")
dev.off()

print("¡Gráfico generado exitosamente!")