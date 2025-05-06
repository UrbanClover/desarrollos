.libPaths(c("/mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/R/R-4.3.0/library", .libPaths()))

# 1. Instalar paquetes --------------------------------------------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.18", ask = FALSE)

BiocManager::install(
  c("AnnotationDbi", "GenomicFeatures", "tximport", "DESeq2", "pheatmap", 
    "ggrepel", "org.Hs.eg.db", "clusterProfiler", "enrichplot"),
  ask = FALSE,
  force = TRUE
)

install.packages("randomForestSRC")
install.packages("tidyverse")

# 2. Cargar bibliotecas -------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicFeatures)
  library(tximport)
  library(DESeq2)
  library(pheatmap)
  library(ggrepel)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(randomForestSRC)
})

# 3. Configurar datos Salmon --------------------------------------------------
muestras <- c("CD1", "CD2", "CD3", "CD4", "Ctrl1", "Ctrl2", "Ctrl3", "Ctrl4")
files <- file.path("resultados", "salmon", muestras, "quant.sf")
names(files) <- muestras
stopifnot(all(file.exists(files)))

# 4. Crear tx2gene ------------------------------------------------------------
txdb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.110.gtf.gz", format = "gtf")
tx2gene <- transcripts(txdb, columns = c("tx_name", "gene_id")) %>% 
  as_tibble() %>%
  transmute(
    transcript_id = as.character(tx_name),
    gene_id = as.character(gene_id)
  ) %>%
  distinct()
write_csv(tx2gene, "tx2gene.csv")

# 5. Importar conteos ---------------------------------------------------------
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                countsFromAbundance = "no", ignoreTxVersion = TRUE)

# 6. Preparar DESeq2 ----------------------------------------------------------
dds <- DESeqDataSetFromTximport(txi, 
                                colData = data.frame(
                                  condition = factor(rep(c("CD", "Ctrl"), each = 4)),
                                  row.names = muestras
                                ),
                                design = ~ condition)

# 7. Filtrado y análisis ------------------------------------------------------
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep,]
dds <- DESeq(dds)

# Calcular resultados DESeq2 (¡Falta clave aquí!)
res <- results(dds, 
              contrast = c("condition", "CD", "Ctrl"), 
              alpha = 0.1,
              lfcThreshold = 1)

vsd <- vst(dds, blind = FALSE)

# 8. Resultados ---------------------------------------------------------------
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene_id") %>%
  mutate(
    significant = ifelse(padj < 0.1 & abs(log2FoldChange) > 1, "Sí", "No"),
    gene_symbol = suppressWarnings(
      mapIds(org.Hs.eg.db,
             keys = gene_id,
             column = "SYMBOL",
             keytype = "ENSEMBL",
             multiVals = function(x) paste(unique(x), collapse = "/"))
    )
  ) %>%
  filter(!is.na(padj))
write_csv(res_df, "Resultados_Completos_Anotados.csv")

# 9. Generar PDF --------------------------------------------------------------
pdf("Informe_Resultados.pdf", width = 12, height = 8)

# PCA Plot (sin cambios)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
print(
  ggplot(pca_data, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(aes(label = name), show.legend = FALSE) +
    ggtitle("Análisis PCA") +
    theme_minimal()
)

# MA Plot mejorado
try({
  ma_data <- res_df %>% arrange(padj)
  print(
    ggplot(ma_data, aes(baseMean, log2FoldChange)) +
      geom_point(aes(color = significant), alpha = 0.5) +
      scale_x_log10() +
      scale_color_manual(values = c("gray70", "#D55E00")) +
      geom_hline(yintercept = 0, color = "grey30") +
      geom_smooth(se = FALSE, color = "navy") +
      ggrepel::geom_text_repel(
        data = subset(ma_data, significant == "Sí" & baseMean > median(baseMean)),
        aes(label = gene_symbol),
        max.overlaps = 100,
        force = 20,
        size = 2.5
      ) +
      labs(
        title = "MA Plot: Expresión vs Cambio Logarítmico",
        x = "Expresión media (log10)",
        y = "Log2 Fold Change"
      ) +
      theme_bw()
  )
})

# Heatmap de genes variables (versión corregida)
tryCatch({
  top_var <- order(rowVars(assay(vsd)), decreasing = TRUE)[1:50]
  annotation_col <- data.frame(
    Condition = colData(vsd)$condition,
    row.names = colnames(vsd)
  )
  
  # Definir colores manualmente
  ann_colors <- list(
    Condition = c(CD = "#1B9E77", Ctrl = "#D95F02")
  )
  
  pheatmap(
    assay(vsd)[top_var, ],
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    main = "50 genes con mayor varianza",
    show_rownames = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
  )
}, error = function(e) {
  message("Error en heatmap: ", e$message)
  plot.new()
  title(main = "Error generando heatmap de variabilidad")
})

# Volcano plot (versión robustecida)
try({
  volcano_data <- res_df %>%
    mutate(
      gene_label = ifelse(significant == "Sí" & !is.na(gene_symbol), 
                          gene_symbol, 
                          "")
    )
  
  print(ggplot(res_df, aes(log2FoldChange, -log10(padj))) +
  geom_point(
    aes(color = significant), 
    alpha = 0.5, 
    size = 1.5,
    shape = 16
  ) +
  scale_color_manual(values = c("gray80", "#D55E00")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey30") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "grey30") +
  ggrepel::geom_text_repel(
    data = subset(res_df, significant == "Sí" & !is.na(gene_symbol)),
    aes(label = gene_symbol),
    max.overlaps = 100,
    force = 20,
    size = 2.8,
    box.padding = 0.15,
    segment.color = "grey50",
    max.time = 10
  ) +
  labs(
    title = "Volcano Plot: CD vs Ctrl",
    subtitle = "Genes significativos (padj < 0.1, |LFC| > 1)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
  )
})


try({
  # Convertir IDs a Entrez
  entrez_ids <- mapIds(org.Hs.eg.db,
                      keys = res_df$gene_id,
                      column = "ENTREZID",
                      keytype = "ENSEMBL",
                      multiVals = "first")
  
  # Enriquecimiento GO
  go_enrich <- enrichGO(
    gene = na.omit(entrez_ids[res_df$significant == "Sí"]),
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
  
  if(nrow(go_enrich) > 0) {
    print(
      dotplot(go_enrich, showCategory=15) +
        ggtitle("Enriquecimiento en términos GO (Biológicos)")
    )
  }
  
  # Enriquecimiento KEGG
  kegg_enrich <- enrichKEGG(
    gene = na.omit(entrez_ids[res_df$significant == "Sí"]),
    organism = 'hsa',
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
  
  if(nrow(kegg_enrich) > 0) {
    print(
      dotplot(kegg_enrich, showCategory=15) +
        ggtitle("Enriquecimiento en rutas KEGG")
    )
  }
})


# Estadísticas resumidas
res_summary <- res_df %>% 
  summarise(
    Genes_Total = n(),
    Genes_Significativos = sum(significant == "Sí"),
    Upregulated = sum(significant == "Sí" & log2FoldChange > 0),
    Downregulated = sum(significant == "Sí" & log2FoldChange < 0),
    Media_Log2FC = mean(log2FoldChange[significant == "Sí"], na.rm = TRUE)
  )

gridExtra::grid.table(t(res_summary), 
                     rows = rownames(t(res_summary)),
                     theme = gridExtra::ttheme_minimal())


dev.off()

# ---------------------------------------------------------------
# 10. Random Forest (Versión Estable)
# ---------------------------------------------------------------
cat("\n10. Análisis con Random Forest\n")

try({
  sig_genes <- res_df %>% 
    filter(significant == "Sí", gene_id %in% rownames(assay(vsd))) %>% 
    pull(gene_id) %>% unique()
  
  if(length(sig_genes) >= 5) {
    expr_matrix <- assay(vsd)[sig_genes, , drop = FALSE]
    rf_data <- data.frame(
      t(expr_matrix),
      Condition = factor(colData(vsd)$condition),
      check.names = FALSE,
      stringsAsFactors = TRUE
    )
    
    # Entrenamiento robusto
    rf_model <- tryCatch({
      rfsrc(Condition ~ ., 
            data = rf_data,
            ntree = 500,
            importance = FALSE,
            seed = 123)
    }, error = function(e) {
      message("Error en RF: ", e$message)
      NULL
    })
    
    if(!is.null(rf_model)) {
      # Cálculo de importancia con verificación
      vimp_result <- tryCatch({
        vimp(rf_model, method = "permute")
      }, error = function(e) {
        message("Error en VIMP: ", e$message)
        NULL
      })
      
      if(!is.null(vimp_result)) {
        vimp_df <- data.frame(
          GeneID = colnames(rf_data)[-ncol(rf_data)],
          Importance = as.numeric(vimp_result$importance),
          Gene_Symbol = suppressWarnings(
            mapIds(org.Hs.eg.db,
                   keys = sub("^X", "", colnames(rf_data)[-ncol(rf_data)]),
                   column = "SYMBOL",
                   keytype = "ENSEMBL",
                   multiVals = function(x) paste(unique(x), collapse = "/"))
          ),
          stringsAsFactors = FALSE,
          row.names = NULL
        ) %>% arrange(desc(Importance))
        
        write_csv(vimp_df, "resultados/importancia_genes.csv")
      }
    }
    
    # ---------------------------------------------------------------
    # 11. Análisis Avanzado (Versión Final)
    # ---------------------------------------------------------------
    cat("\n11. Análisis Avanzado\n")
    
    # Matriz de proximidad con nombres
    if(!is.null(rf_model$proximity)) {
      pheatmap(rf_model$proximity,
               main = "Matriz de Proximidad",
               color = colorRampPalette(c("#F0F9E8", "#0868AC"))(100),
               show_rownames = TRUE,
               show_colnames = TRUE)
    }
    
    # Error OOB con datos validados
    if(!is.null(rf_model$err.rate)) {
      error_df <- data.frame(
        Árboles = seq_len(nrow(rf_model$err.rate)),
        Error_OOB = rf_model$err.rate[,1]
      ) %>% na.omit()
      
      if(nrow(error_df) > 1) {
        print(
          ggplot(error_df, aes(Árboles, Error_OOB)) +
            geom_line(color = "#D73027", linewidth = 1) +
            labs(title = "Error OOB") +
            theme_bw()
        )
      }
    }
    
    # Supervivencia con datos robustos
    surv_data <- data.frame(
      Tiempo = runif(nrow(rf_data), 30, 365),
      Estado = sample(0:1, nrow(rf_data), replace = TRUE),  # replace=TRUE crítico
      rf_data[, -ncol(rf_data)],
      check.names = FALSE
    )
    
    surv_model <- tryCatch({
      rfsrc(Surv(Tiempo, Estado) ~ ., 
            data = surv_data,
            ntree = 500)
    }, error = function(e) {
      message("Error en supervivencia: ", e$message)
      NULL
    })
  }
})