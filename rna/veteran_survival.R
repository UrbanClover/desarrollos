.libPaths(c("/mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/R/R-4.4.3-build/lib64", .libPaths()))

# ------------------------------------------------------------
# PASO 0: Instalar y cargar paquetes
# ------------------------------------------------------------
if (!require("survival")) install.packages("survival")
if (!require("randomForestSRC")) install.packages("randomForestSRC")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("caTools")) install.packages("caTools")
if (!require("pec")) install.packages("pec")          # Para Brier Score
if (!require("timeROC")) install.packages("timeROC")  # Para ROC temporal
if (!require("fastshap")) install.packages("fastshap")# Para SHAP values
if (!require("survminer")) install.packages("survminer")

library(survival)
library(randomForestSRC)
library(ggplot2)
library(caTools)
library(pec)
library(timeROC)
library(fastshap)
library(survminer)

# ------------------------------------------------------------
# PASO 1: Cargar y preparar los datos (veteran)
# ------------------------------------------------------------
data(veteran)

# Convertir variables categóricas a factor
veteran$celltype <- as.factor(veteran$celltype)
veteran$trt <- as.factor(veteran$trt)

# Eliminar filas con NA
veteran_clean <- na.omit(veteran)

# ------------------------------------------------------------
# PASO 2: Dividir en entrenamiento (80%) y prueba (20%)
# ------------------------------------------------------------
set.seed(123)
split <- sample.split(veteran_clean$status, SplitRatio = 0.8)
train_data <- subset(veteran_clean, split == TRUE)
test_data <- subset(veteran_clean, split == FALSE)

cat("Tamaño entrenamiento:", nrow(train_data), "\n")
cat("Tamaño prueba:", nrow(test_data), "\n")

# ------------------------------------------------------------
# PASO 3: Entrenar el modelo
# ------------------------------------------------------------
set.seed(123)
rsf_model <- rfsrc(
  Surv(time, status) ~ .,
  data = train_data,
  ntree = 2000,
  nodesize = 10,
  importance = TRUE,
  samptype = "swor"
)

# ------------------------------------------------------------
# PASO 4: Evaluación del modelo
# ------------------------------------------------------------
# 1. Error OOB y VIMP
cat("Error OOB (entrenamiento):", rsf_model$err.rate[rsf_model$ntree], "\n")

vimp <- rsf_model$importance
vimp_df <- data.frame(Variable = names(vimp), Importance = vimp)

ggplot(vimp_df, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  labs(title = "Importancia de Variables (VIMP)", x = "")

# 2. Predicciones en test
pred_test <- predict(rsf_model, newdata = test_data)

# ------------------------------------------------------------
# PASO 5: Métricas avanzadas
# ------------------------------------------------------------
# 1. Brier Score

pred_vars <- setdiff(names(test_data), c("time","status"))

fmla <- as.formula(paste("Surv(time, status) ~", paste(pred_vars, collapse = "+")))

brier <- pec(
  list(RSF = rsf_model),
  formula = fmla,
  data = test_data,
  times = quantile(test_data$time, probs = c(0.25, 0.5, 0.75))
)
print(brier)

# 2. ROC Temporal
roc_data <- timeROC(
  T = test_data$time,
  delta = test_data$status,
  marker = pred_test$predicted,
  cause = 1,
  times = c(100, 200, 300),
  ROC = TRUE
)
plot(roc_data, time = 200, title = "ROC a 200 días")

# 3. Comparación con Cox
cox_model <- coxph(Surv(time, status) ~ ., data = train_data)
pred_cox <- predict(cox_model, newdata = test_data, type = "risk")
cindex_cox <- concordance(Surv(time, status) ~ pred_cox, data = test_data, reverse = TRUE)$concordance
cat("C-index (Cox):", cindex_cox, "\n")

# 4. SHAP Values (Opcional: Requiere tiempo de cómputo)
pred_wrapper <- function(model, newdata) predict(model, newdata)$predicted
shap_values <- explain(
  rsf_model, 
  X = train_data[,!names(train_data) %in% c("time","status")], 
  pred_wrapper = pred_wrapper,
  nsim=10)  # nsim reducido para prueba
plot(shap_values)

# 5. Clustering de riesgo
risk_groups <- kmeans(pred_test$predicted, centers = 2)
test_data$risk_group <- as.factor(risk_groups$cluster)
surv_fit <- survfit(Surv(time, status) ~ risk_group, data = test_data)
ggsurvplot(surv_fit, data = test_data, pval = TRUE, risk.table = TRUE)

# ------------------------------------------------------------
# PASO 6: Guardar modelo
# ------------------------------------------------------------
saveRDS(rsf_model, "modelo_entrenado_80_20.RDS")