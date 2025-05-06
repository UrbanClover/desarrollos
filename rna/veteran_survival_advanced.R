# ------------------------------------------------------------
# CONFIGURACIÓN INICIAL
# ------------------------------------------------------------
.libPaths(c("/mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/R/R-4.4.3-build/lib64", .libPaths()))

# Instalar y cargar paquetes
required_packages <- c("survival", "randomForestSRC", "ggplot2", "caTools", "dplyr",
                      "Hmisc", 
                      "pec", "timeROC", "fastshap", "survminer", "coxphf", "report")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# ------------------------------------------------------------
# PASO 1: Carga y preparación de datos
# ------------------------------------------------------------
data(veteran)
veteran <- veteran %>%
  dplyr::mutate(
    celltype = factor(celltype),
    trt = factor(trt),
    prior = factor(prior)
  ) %>%
  na.omit()

set.seed(123)
split <- sample.split(veteran$status, SplitRatio = 0.8)
train_data <- subset(veteran, split)
test_data <- subset(veteran, !split)

# ------------------------------------------------------------
# PASO 2: Entrenamiento del modelo con optimización
# ------------------------------------------------------------
# Búsqueda de hiperparámetros
param_grid <- expand.grid(
  nodesize = c(5, 10),
  mtry = c(3, 4),
  ntree = c(1000, 2000)
)

best_error <- Inf
best_model <- NULL

for(i in 1:nrow(param_grid)) {
  cat("\nProbando combinación", i, "/", nrow(param_grid), "...\n")
  current_model <- rfsrc(
    Surv(time, status) ~ .,
    data = train_data,
    ntree = param_grid$ntree[i],
    mtry = param_grid$mtry[i],
    nodesize = param_grid$nodesize[i],
    importance = TRUE
  )
  
  if(current_model$err.rate[current_model$ntree] < best_error) {
    best_model <- current_model
    best_error <- current_model$err.rate[current_model$ntree]
    best_params <- param_grid[i,]
  }
}

rsf_model <- best_model
cat("\nMejores parámetros encontrados:\n")
print(best_params)

# ------------------------------------------------------------
# PASO 3: Diagnóstico del modelo
# ------------------------------------------------------------
# 1. Convergencia del modelo
png("convergence_plot.png", width=800, height=600)
plot(rsf_model$err.rate, type = "l", lwd = 2,
     main = "Error OOB durante el entrenamiento",
     xlab = "Número de árboles", ylab = "Error OOB",
     col = "darkred")
dev.off()

# 2. Importancia de variables
vimp_df <- data.frame(
  Variable = names(rsf_model$importance),
  Importance = rsf_model$importance
) %>% arrange(desc(Importance))

ggplot(vimp_df, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Importancia de Variables (VIMP)", 
       subtitle = paste("Error OOB final:", round(best_error, 4)),
       x = "") +
  theme_minimal()

ggsave("vimp_plot.png", width = 10, height = 6)

# ------------------------------------------------------------
# PASO 4: Evaluación del rendimiento
# ------------------------------------------------------------
# Predicciones en test
pred_test <- predict(rsf_model, newdata = test_data)

# 1. Métricas comparativas
metrics <- list(
  RSF = list(
    Cindex = Hmisc::rcorr.cens(pred_test$predicted, Surv(test_data$time, test_data$status))["C Index"],
    Brier = pec(list(RSF = rsf_model), 
               formula = Surv(time, status) ~ ., 
               data = test_data)$AppErr$RSF
  ),
  Cox = {
    cox_model <- coxphf(Surv(time, status) ~ ., data = train_data)
    pred_cox <- predict(cox_model, newdata = test_data, type = "risk")
    list(
      Cindex = concordance(Surv(time, status) ~ pred_cox, 
                          data = test_data, reverse = TRUE)$concordance,
      Brier = pec(list(Cox = cox_model), 
                 formula = Surv(time, status) ~ ., 
                 data = test_data)$AppErr$Cox
    )
  }
)

# Tabla comparativa
metrics_table <- data.frame(
  Modelo = c("RSF", "Cox"),
  Cindex = c(metrics$RSF$Cindex, metrics$Cox$Cindex),
  BrierScore = c(tail(metrics$RSF$Brier, 1), tail(metrics$Cox$Brier, 1))
)

print(metrics_table)

# 2. Curvas de calibración
cal_plot <- calibration_plot(
  predicted = pred_test$predicted,
  time = test_data$time,
  status = test_data$status,
  groups = 4
) + 
  labs(title = "Calibración del Modelo RSF",
       subtitle = "Comparación entre riesgo predicho y observado")

ggsave("calibration_plot.png", cal_plot, width = 10, height = 6)

# ------------------------------------------------------------
# PASO 5: Análisis de supervivencia avanzado
# ------------------------------------------------------------
# 1. Grupos de riesgo
risk_groups <- cut(pred_test$predicted, 
                  breaks = quantile(pred_test$predicted, probs = c(0, 0.33, 0.66, 1)),
                  labels = c("Bajo", "Medio", "Alto"))

surv_fit <- survfit(Surv(time, status) ~ risk_groups, data = test_data)

surv_plot <- ggsurvplot(
  surv_fit,
  data = test_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = "abs_pct",
  palette = "jco",
  title = "Supervivencia por Grupos de Riesgo",
  legend.labs = levels(risk_groups)
)

ggsave("survival_plot.png", print(surv_plot), width = 12, height = 8)

# 2. Análisis de subgrupos
celltype_plot <- ggsurvplot(
  survfit(Surv(time, status) ~ celltype, data = test_data),
  data = test_data,
  pval = TRUE,
  conf.int = TRUE,
  palette = "Dark2",
  title = "Supervivencia por Tipo Celular"
)

ggsave("celltype_plot.png", print(celltype_plot), width = 12, height = 8)

# ------------------------------------------------------------
# PASO 6: Interpretación del modelo
# ------------------------------------------------------------
# 1. Análisis SHAP (versión optimizada)
set.seed(123)
shap_values <- explain(
  rsf_model, 
  X = train_data[,!names(train_data) %in% c("time","status")], 
  nsim = 20,  # Reducido para eficiencia
  .progress = "text",
  pred_wrapper = function(model, newdata) predict(model, newdata)$predicted
)

# Gráficos SHAP
shap_summary <- plot(shap_values, type = "importance", num_features = 15) +
  labs(title = "Importancia SHAP",
       subtitle = "Contribución promedio absoluta de cada variable")

shap_dependence <- plot(shap_values, type = "dependence", 
                       feature = "karno", X = train_data) +
  labs(title = "Dependencia SHAP para Karnofsky Score",
       y = "Valor SHAP")

ggsave("shap_importance.png", shap_summary, width = 10, height = 6)
ggsave("shap_dependence.png", shap_dependence, width = 10, height = 6)

# ------------------------------------------------------------
# PASO 7: Validación y reporte final
# ------------------------------------------------------------
# 1. Validación bootstrap
boot_results <- validate(
  rsf_model,
  method = "boot",
  B = 50,  # Reducido para eficiencia
  dxy = TRUE,
  metrics = c("Cindex", "Brier")
)

# 2. Reporte HTML automatizado
report_content <- list(
  Model = rsf_model,
  Metrics = metrics,
  SHAP = shap_values,
  Validation = boot_results
)

html_report <- report(report_content)
export(html_report, "survival_analysis_report.html")

# 3. Guardar todos los resultados
saveRDS(list(
  model = rsf_model,
  metrics = metrics,
  shap_values = shap_values,
  validation = boot_results
), "full_analysis_results.RDS")

# Mensaje final
cat("\nAnálisis completado con éxito!\n",
    "Resultados guardados en:\n",
    "- survival_analysis_report.html\n",
    "- full_analysis_results.RDS\n",
    "- Varios archivos .png con gráficos\n")