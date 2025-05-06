# ------------------------------------------------------------
# CONFIGURACIÓN INICIAL
# ------------------------------------------------------------
.libPaths(c("/mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/R/R-4.4.3-build/lib64", .libPaths()))

# Instalar paquetes (si no están instalados)
install.packages(c("randomForest", "caret", "ggplot2", "pROC", "dplyr"),
                 lib = "/mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/R/R-4.4.3-build/lib64",
                 repos = "https://cloud.r-project.org",
                 dependencies = TRUE,
                 verbose = TRUE)

# Cargar librerías
library(randomForest)  # Modelo Random Forest
library(caret)         # Preprocesamiento y métricas
library(ggplot2)       # Visualización
library(pROC)          # Curva ROC
library(dplyr)         # Manipulación de datos

# ------------------------------------------------------------

# Cargar el dataset (ajusta la ruta a tu archivo)
datos <- read.csv("diabetes.csv")

# Ver estructura (asegúrate que la variable objetivo es 'Outcome')
str(datos)

# Convertir 'Outcome' a factor y renombrar niveles para claridad
datos$Outcome <- as.factor(datos$Outcome)
levels(datos$Outcome) <- c("No_Diabetico", "Diabetico")

# Verificar valores faltantes o ceros inválidos (ej: glucosa = 0)
# En este dataset, los ceros en variables como 'Glucose' son errores. Los reemplazamos por NA:
datos <- datos %>%
  mutate(
    Glucose = ifelse(Glucose == 0, NA, Glucose),
    BloodPressure = ifelse(BloodPressure == 0, NA, BloodPressure),
    SkinThickness = ifelse(SkinThickness == 0, NA, SkinThickness),
    Insulin = ifelse(Insulin == 0, NA, Insulin),
    BMI = ifelse(BMI == 0, NA, BMI)
  )

# Imputar valores faltantes con la mediana (usando caret)
preprocesador <- preProcess(datos, method = c("medianImpute"))
datos <- predict(preprocesador, datos)

# Verificar que no hay NA
sum(is.na(datos))

# ------------------------------------------------------------

# DIVIDIR EN TRAIN Y TEST

set.seed(123)  # Para reproducibilidad
indices <- createDataPartition(datos$Outcome, p = 0.8, list = FALSE)
entrenamiento <- datos[indices, ]
prueba <- datos[-indices, ]

# ------------------------------------------------------------

# MODELO RANDOM FOREST

# Configurar validación cruzada (10-fold)
ctrl <- trainControl(
  method = "cv",          # Validación cruzada simple
  number = 10,            # 10 particiones
  classProbs = TRUE,      # Necesario para la curva ROC
  summaryFunction = twoClassSummary  # Métricas para clases
)

# Entrenar el modelo
modelo <- train(
  Outcome ~ ., 
  data = entrenamiento,
  method = "rf",
  trControl = ctrl,
  tuneLength = 5,         # Ajusta automáticamente 'mtry'
  metric = "ROC"          # Optimizar por AUC-ROC
)

# Ver resultados del modelo
print(modelo)
plot(modelo)  # Ver rendimiento según 'mtry'

# ------------------------------------------------------------

# EVALUACIÓN DEL MODELO

# Predecir en datos de prueba
predicciones <- predict(modelo, newdata = prueba)

# Matriz de confusión (detalle de aciertos/errores)
confusionMatrix(predicciones, prueba$Outcome)

# Curva ROC y AUC
probabilidades <- predict(modelo, prueba, type = "prob")[, "Diabetico"]


# PASO 5: GUARDAR GRÁFICOS EN UN PDF CON NOMBRE PERSONALIZADO
pdf("Resultados_RandomForest_Diabetes.pdf")  # Abrir dispositivo PDF

# Gráfico 1: Curva ROC
roc_curve <- roc(prueba$Outcome, probabilidades, levels = c("No_Diabetico", "Diabetico"))
plot(roc_curve, main = "Curva ROC - Diabetes")

# Gráfico 2: Importancia de variables
importancia <- varImp(modelo)
ggplot(importancia, top = 8) + 
  geom_bar(stat = "identity", fill = "steelblue") +
  ggtitle("Importancia de variables en el modelo") +
  theme_minimal()

dev.off()  # Cerrar y guardar el PDF

# ------------------------------------------------------------
# MODIFICACIÓN: AÑADIR ANÁLISIS DE SUPERVIVENCIA
# ------------------------------------------------------------

# Instalar paquetes adicionales si son necesarios
install.packages(c("survival", "survminer"),
                 lib = "/mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/R/R-4.4.3-build/lib64",
                 repos = "https://cloud.r-project.org",
                 dependencies = TRUE,
                 verbose = TRUE)

# Cargar nuevas librerías
library(survival)     # Análisis de supervivencia
library(survminer)    # Visualización de supervivencia

# ------------------------------------------------------------
# SIMULAR DATOS DE SUPERVIVENCIA (FICTICIOS)
# ------------------------------------------------------------

set.seed(123)
# Crear variables ficticias para demostración:
datos <- datos %>%
  mutate(
    # Tiempo hasta evento (en años): valores entre 1-15
    survival_time = runif(n(), 1, 15),
    
    # Evento (0 = censurado, 1 = muerte): 
    # Relacionado ficticiamente con el Outcome
    survival_event = ifelse(
      Outcome == "Diabetico", 
      rbinom(n(), 1, 0.7),  # Mayor probabilidad de evento en diabéticos
      rbinom(n(), 1, 0.3)    # Menor probabilidad en no diabéticos
    )
  )

# Ver estructura de los nuevos datos
str(datos[c("survival_time", "survival_event")])

# ------------------------------------------------------------
# ANÁLISIS DE SUPERVIVENCIA
# ------------------------------------------------------------

# 1. Curva de Kaplan-Meier (supervivencia general)
km_fit <- survfit(Surv(survival_time, survival_event) ~ 1, data = datos)
print(km_fit)

# 2. Curvas de supervivencia por grupo (Diabéticos vs No Diabéticos)
km_group_fit <- survfit(Surv(survival_time, survival_event) ~ Outcome, data = datos)

# Gráfico comparativo
km_plot <- ggsurvplot(
  km_group_fit,
  data = datos,
  pval = TRUE,          # Muestra p-valor de log-rank test
  conf.int = TRUE,      # Intervalos de confianza
  palette = c("#E7B800", "#2E9FDF"),  # Colores personalizados
  xlab = "Tiempo (años)", 
  ylab = "Probabilidad de Supervivencia",
  legend.title = "Grupo",
  risk.table = TRUE     # Tabla de riesgo añadida
)

# 3. Modelo de Cox (efecto de variables en el riesgo)
cox_model <- coxph(
  Surv(survival_time, survival_event) ~ Age + BMI + Glucose + Outcome,
  data = datos
)

# Resumen del modelo
summary(cox_model)

# ------------------------------------------------------------
# ACTUALIZAR EL PDF CON GRÁFICOS DE SUPERVIVENCIA
# ------------------------------------------------------------

pdf("Resultados_RandomForest_Diabetes_CON_Supervivencia.pdf")

# Gráficos originales
plot(roc_curve, main = "Curva ROC - Diabetes")
ggplot(importancia, top = 8) + 
  geom_bar(stat = "identity", fill = "steelblue") +
  ggtitle("Importancia de variables en el modelo") +
  theme_minimal()

# Nuevos gráficos de supervivencia
print(km_plot, newpage = FALSE)  # Curvas de supervivencia por grupo

# Gráfico del modelo de Cox
ggforest(cox_model, data = datos)  # Efectos de las variables en el riesgo

dev.off()

# ------------------------------------------------------------