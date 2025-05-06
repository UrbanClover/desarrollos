.libPaths(c("/mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/R/R-4.3.0/library", .libPaths()))

library(qqman)
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)

# Leer datos y filtrar solo filas con TEST == "ADD"
gwas_data <- read.table(args[1], header = TRUE, sep = "\t", comment.char = "")
gwas_data <- subset(gwas_data, TEST == "ADD")  # ¡Filtrar solo ADD!

# Convertir P a numérico y eliminar NA
gwas_data$P <- as.numeric(gwas_data$P)
gwas_data <- na.omit(gwas_data)  # Eliminar filas con P=NA

# Generar QQ plot
png(args[2])
qq(gwas_data$P, main = "QQ Plot")
dev.off()