# ------------------------------------------------------------
# CONFIGURACIÓN INICIAL
# ------------------------------------------------------------
.libPaths(c("/mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/R/R-4.4.3-build/lib64", .libPaths()))

# Instalar paquetes (si no están instalados)
install.packages(c("randomForest", "caret", "ggplot2", "pROC", "dplyr", "xgboost", "reshape2", "kernlab", "gridExtra"),
                 lib = "/mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/R/R-4.4.3-build/lib64",
                 repos = "https://cloud.r-project.org",
                 dependencies = TRUE,
                 verbose = TRUE)

# ------------------------------------------------------------
# CONFIGURACIÓN INICIAL (SIMPLIFICADA)
# ------------------------------------------------------------
# Cargar librerías sin instalación (asumiendo que ya están instaladas)
library(randomForest)
library(caret)
library(pROC)
library(dplyr)
library(xgboost)
library(ggplot2)
library(gridExtra)

nombre_archivo <- "Analisis_Diabetes_Modelos.pdf"

# ------------------------------------------------------------
# PASO 2: CARGAR Y PREPROCESAR DATOS (SIMPLIFICADO)
# -----------------------------------
datos <- read.csv("diabetes.csv") %>% 
  mutate(
    Outcome = factor(Outcome, levels = c(0, 1), labels = c("No_Diabetico", "Diabetico"))
  ) %>% 
  na.omit()  # Eliminar filas con NA en lugar de imputar

# Dividir en entrenamiento (80%) y prueba (20%)
set.seed(123)
indices <- createDataPartition(datos$Outcome, p = 0.8, list = FALSE)
entrenamiento <- datos[indices, ]
prueba <- datos[-indices, ]

# ------------------------------------------------------------
# PASO 3: ENTRENAR MODELOS ESENCIALES (EVITAR SVM)
# ------------------------------------------
ctrl <- trainControl(
  method = "cv",  # Validación cruzada simple
  number = 5,     # 5 folds para mayor velocidad
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  verboseIter = TRUE
)

# Modelo 1: Random Forest
modelo_rf <- train(
  Outcome ~ .,
  data = entrenamiento,
  method = "rf",
  trControl = ctrl,
  tuneLength = 5,
  metric = "ROC",
  importance = TRUE,
  probability = TRUE
)

# Modelo 2: XGBoost
modelo_xgb <- train(
  Outcome ~ .,
  data = entrenamiento,
  method = "xgbTree",
  trControl = ctrl,
  tuneLength = 5,
  metric = "ROC",
  objective = "binary:logistic",
  eval_metric = "auc",
  verbose = 0
)

# ------------------------------------------------------------
# PASO 4: EVALUACIÓN DIRECTA (VERSIÓN CORREGIDA)
# ---------------------------------------------
# Función para obtener ROC y métricas
evaluar_modelo <- function(modelo) {
  prob <- predict(modelo, prueba, type = "prob")[, "Diabetico"]
  pred_clases <- predict(modelo, prueba)
  
  # ROC y AUC
  roc_obj <- roc(prueba$Outcome, prob)
  auc_val <- auc(roc_obj)
  
  # Matriz de confusión
  cm <- confusionMatrix(pred_clases, prueba$Outcome, positive = "Diabetico")
  
  return(list(
    ROC = roc_obj,
    AUC = auc_val,
    Sensibilidad = cm$byClass["Sensitivity"],
    Especificidad = cm$byClass["Specificity"],
    MatrizConfusion = cm
  ))
}

# Generar objetos para gráficos
resultados_rf <- evaluar_modelo(modelo_rf)
resultados_xgb <- evaluar_modelo(modelo_xgb)

# Tabla comparativa
tabla_comparativa <- data.frame(
  Modelo = c("Random Forest", "XGBoost"),
  Sensibilidad = c(resultados_rf$Sensibilidad, resultados_xgb$Sensibilidad),
  Especificidad = c(resultados_rf$Especificidad, resultados_xgb$Especificidad),
  AUC = c(resultados_rf$AUC, resultados_xgb$AUC)
)
print(tabla_comparativa)

# ------------------------------------------------------------
# FUNCIONES PARA GRÁFICOS (DEFINIR ANTES DE USARLAS)
# ------------------------------------------------------------
plot_confusion_matrix <- function(cm, titulo) {
  df <- as.data.frame(cm$table)
  ggplot(df, aes(Prediction, Reference, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq), color = "white", size = 5) +
    scale_fill_gradient(low = "blue", high = "red") +
    labs(title = titulo, x = "Predicción", y = "Real") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
}

generar_curva_aprendizaje <- function(modelo, datos, metrica = "ROC") {
  curva <- caret::learning_curve_dat(
    datos,
    outcome = "Outcome",
    metric = metrica,
    test_prop = 0.2,
    verbose = FALSE
  )
  ggplot(curva, aes(x = Training_Size, y = Metric, color = Data)) +
    geom_smooth(method = "loess", span = 0.8) +
    labs(title = paste("Curva de Aprendizaje -", modelo$modelInfo$label),
         x = "Tamaño del Conjunto de Entrenamiento",
         y = metrica) +
    theme_minimal()
}

# ------------------------------------------------------------
# GENERACIÓN DE OBJETOS PARA GRÁFICOS (DESPUÉS DE DEFINIR FUNCIONES)
# ------------------------------------------------------------
curva_rf <- generar_curva_aprendizaje(modelo_rf, entrenamiento)
curva_xgb <- generar_curva_aprendizaje(modelo_xgb, entrenamiento)
tabla_melt <- reshape2::melt(tabla_comparativa, id.vars = "Modelo")

# ------------------------------------------------------------
# GENERACIÓN DE GRÁFICOS (VERSIÓN FINAL FUNCIONAL)
# ------------------------------------------------------------
pdf(nombre_archivo, width = 14, height = 10)

# 1. Curvas ROC Comparativas
plot(resultados_rf$ROC, col = "red", main = "Curvas ROC Comparativas")
lines(resultados_xgb$ROC, col = "blue")
legend("bottomright", 
       legend = c(paste("RF (AUC =", round(resultados_rf$AUC, 2)),
                  paste("XGB (AUC =", round(resultados_xgb$AUC, 2))),
       col = c("red", "blue"), lwd = 2)

# 2. Matrices de Confusión
grid.arrange(
  plot_confusion_matrix(resultados_rf$MatrizConfusion, "Random Forest") + theme(legend.position = "none"),
  plot_confusion_matrix(resultados_xgb$MatrizConfusion, "XGBoost") + theme(legend.position = "none"),
  ncol = 2
)

# 3. Importancia de Variables
grid.arrange(
  plot_importancia_variables(modelo_rf, "Random Forest"),
  plot_importancia_variables(modelo_xgb, "XGBoost"),
  ncol = 2
)

# 4. Distribución de Probabilidades
grid.arrange(
  plot_distribucion_prob(predict(modelo_rf, prueba, type = "prob")[,2], prueba$Outcome),
  plot_distribucion_prob(predict(modelo_xgb, prueba, type = "prob")[,2], prueba$Outcome),
  ncol = 2
)

# 5. Heatmap de Correlación
print(plot_heatmap_correlacion(datos))

# 6. Curvas de Aprendizaje
grid.arrange(
  curva_rf + theme(legend.position = "none"),
  curva_xgb + theme(legend.position = "none"),
  ncol = 2
)

# 7. Comparación de Métricas
print(
  ggplot(tabla_melt, aes(x = Modelo, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = "Comparación de Métricas por Modelo", y = "Valor") +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal() +
    theme(text = element_text(size = 12))
)

dev.off()