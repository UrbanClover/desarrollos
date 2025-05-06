# ------------------------------------------------------------
# CONFIGURACIÓN INICIAL
# ------------------------------------------------------------
.libPaths(c("/mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/R/R-4.4.3-build/lib64", .libPaths()))

# Instalar paquetes (si no están instalados)
install.packages(c("randomForest", "caret", "ggplot2", "pROC", "dplyr", "xgboost", "reshape2", "kernlab"),
                 lib = "/mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/R/R-4.4.3-build/lib64",
                 repos = "https://cloud.r-project.org",
                 dependencies = c("Depends", "Imports", "LinkingTo"),
                 verbose = TRUE)

# Cargar librerías
library(randomForest)  # Random Forest
library(caret)         # Preprocesamiento y evaluación
library(pROC)          # Curvas ROC
library(dplyr)         # Manipulación de datos
library(xgboost)       # XGBoost
library(ggplot2)       # Gráficos
library(reshape2)      # Formateo de datos para gráficos
library(kernlab)

# ------------------------------------------------------------
# PASO 2: CARGAR Y PREPROCESAR DATOS
# -----------------------------------
# Cargar datos (ajusta la ruta a tu archivo)
datos <- read.csv("diabetes.csv")

# Convertir 'Outcome' a factor con etiquetas claras
datos$Outcome <- as.factor(datos$Outcome)
levels(datos$Outcome) <- c("No_Diabetico", "Diabetico")

# Manejar ceros como NA (variables clínicas no pueden ser 0)
datos <- datos %>%
  mutate(
    Glucose = ifelse(Glucose == 0, NA, Glucose),
    BloodPressure = ifelse(BloodPressure == 0, NA, BloodPressure),
    SkinThickness = ifelse(SkinThickness == 0, NA, SkinThickness),
    Insulin = ifelse(Insulin == 0, NA, Insulin),
    BMI = ifelse(BMI == 0, NA, BMI)
  )

# Imputar NA con mediana y escalar (para SVM)
preprocesador <- preProcess(datos, method = c("medianImpute", "center", "scale"))
datos <- predict(preprocesador, datos)

# Dividir en entrenamiento (80%) y prueba (20%)
set.seed(123)
indices <- createDataPartition(datos$Outcome, p = 0.8, list = FALSE)
entrenamiento <- datos[indices, ]
prueba <- datos[-indices, ]


# PASO 3: CONFIGURAR ENTRENAMIENTO Y MODELOS
# ------------------------------------------
# Control de validación cruzada (10-fold)
ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  verboseIter = TRUE,
)

# Entrenar los tres modelos
# --- Modelo 1: Random Forest ---
modelo_rf <- train(
  Outcome ~ .,
  data = entrenamiento,
  method = "rf",
  trControl = ctrl,
  tuneLength = 5,
  metric = "ROC"
)

# --- Modelo 2: SVM Lineal ---
modelo_svm <- train(
  Outcome ~ .,
  data = entrenamiento,
  method = "svmRadial",
  trControl = ctrl,
  tuneLength = 5,
  metric = "ROC",
  prob.model = TRUE
)

library(kernlab)

# --- Modelo 3: XGBoost ---
modelo_xgb <- train(
  Outcome ~ .,
  data = entrenamiento,
  method = "xgbTree",
  trControl = ctrl,
  tuneLength = 5,
  metric = "ROC",
  objective = "binary:logistic"
)

# ------------------------------------------------------------
# PASO 3.1: CREAR OBJETO 'resultados' PARA COMPARACIÓN
# ------------------------------------------------------------
resultados <- resamples(list(
  RF = modelo_rf,
  SVM = modelo_svm,
  XGB = modelo_xgb
))

# PASO 4: EVALUAR MODELOS Y RECOLECTAR MÉTRICAS
# ---------------------------------------------
# Función para extraer métricas
evaluar_modelo <- function(modelo, datos_prueba) {
  # Predecir clases y probabilidades
  pred_clases <- predict(modelo, datos_prueba)
  pred_probs <- predict(modelo, datos_prueba, type = "prob")[, "Diabetico"]
  
  # Matriz de confusión
  cm <- confusionMatrix(pred_clases, datos_prueba$Outcome, positive = "Diabetico")
  
  # Curva ROC y AUC
  roc_obj <- roc(response = datos_prueba$Outcome, predictor = pred_probs, levels = c("No_Diabetico", "Diabetico"))
  auc_val <- auc(roc_obj)
  
  return(list(
    Sensibilidad = cm$byClass["Sensitivity"],
    Especificidad = cm$byClass["Specificity"],
    AUC = auc_val,
    ROC = roc_obj
  ))
}

# Evaluar cada modelo
resultados_rf <- evaluar_modelo(modelo_rf, prueba)
resultados_svm <- evaluar_modelo(modelo_svm, prueba)
resultados_xgb <- evaluar_modelo(modelo_xgb, prueba)

# Crear tabla comparativa
tabla_comparativa <- data.frame(
  Modelo = c("Random Forest", "SVM Lineal", "XGBoost"),
  Sensibilidad = c(resultados_rf$Sensibilidad, resultados_svm$Sensibilidad, resultados_xgb$Sensibilidad),
  Especificidad = c(resultados_rf$Especificidad, resultados_svm$Especificidad, resultados_xgb$Especificidad),
  AUC = c(resultados_rf$AUC, resultados_svm$AUC, resultados_xgb$AUC)
)

print(tabla_comparativa)

# ------------------------------------------------------------
# FUNCIONES AUXILIARES PARA GRÁFICOS (ACTUALIZADO)
# ------------------------------------------------------------
# Función para matriz de confusión
plot_confusion_matrix <- function(cm, titulo) {
  df <- as.data.frame(cm$table)
  ggplot(df, aes(Prediction, Reference, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq), color = "white") +
    scale_fill_gradient(low = "blue", high = "red") +
    labs(title = titulo, x = "Predicción", y = "Real") +
    theme_minimal()
}

# Función para curva de aprendizaje
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

# Generar curvas
curva_rf <- generar_curva_aprendizaje(modelo_rf, entrenamiento)
curva_svm <- generar_curva_aprendizaje(modelo_svm, entrenamiento)
curva_xgb <- generar_curva_aprendizaje(modelo_xgb, entrenamiento)

# Generar matriz de confusión para RF
cm_rf <- confusionMatrix(predict(modelo_rf, prueba), prueba$Outcome)


# PASO 5: GRÁFICOS COMPARATIVOS
pdf("Comparacion_Modelos_Diabetes.pdf", width = 12, height = 8)

# Gráfico 1: Curvas ROC
plot(resultados_rf$ROC, col = "red", main = "Curvas ROC Comparativas")
lines(resultados_svm$ROC, col = "blue")
lines(resultados_xgb$ROC, col = "green")
legend("bottomright", 
       legend = c(paste("RF (AUC =", round(resultados_rf$AUC, 2)),
                  paste("SVM (AUC =", round(resultados_svm$AUC, 2)),
                  paste("XGB (AUC =", round(resultados_xgb$AUC, 2))),
       col = c("red", "blue", "green"), lwd = 2)

# Gráfico 2: Sensibilidad y Especificidad
tabla_melt <- melt(tabla_comparativa, id.vars = "Modelo", measure.vars = c("Sensibilidad", "Especificidad"))
print(
  ggplot(tabla_melt, aes(x = Modelo, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = "Comparación de Sensibilidad y Especificidad", y = "Valor", fill = "Métrica") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set1")
)

# Gráficos nuevos
print(bwplot(resultados, metric = "ROC"))  # Boxplot de AUC
print(learning_curve_xgb)                  # Curva de aprendizaje de XGBoost
print(plot_confusion_matrix(cm_rf, "Matriz de Confusión - Random Forest"))

dev.off()

# PASO 6: SELECCIÓN DEL MEJOR MODELO
# ------------------------------------------------------------
# Obtener métricas promedio
metricas_promedio <- summary(resultados)

# Seleccionar modelo con mayor AUC
mejor_modelo <- names(which.max(metricas_promedio$statistics$ROC[, "Mean"]))
cat("El mejor modelo es:", mejor_modelo)

# Guardar mejor modelo
if (mejor_modelo == "RF") {
  saveRDS(modelo_rf, "mejor_modelo.rds")
} else if (mejor_modelo == "XGB") {
  saveRDS(modelo_xgb, "mejor_modelo.rds")
}