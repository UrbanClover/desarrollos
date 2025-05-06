# --------------------------------------------
# PASO 1: Instalación e carga de paquetes
# --------------------------------------------
if (!require("randomForestSRC")) install.packages("randomForestSRC")
if (!require("survival")) install.packages("survival")
if (!require("survminer")) install.packages("survminer")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")  # Novo paquete para mutate

library(randomForestSRC)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

# --------------------------------------------
# PASO 2: Carga de datos e preparación (CORRECCIÓN)
# --------------------------------------------
data <- read.delim("simulated_patients.tsv", sep = "\t", stringsAsFactors = TRUE)

# Converter variables específicas e character a factor/numeric
data <- data %>%
  mutate(
    sexo = as.factor(sexo),
    mutacion_LTBP4_IDH1 = as.factor(mutacion_LTBP4_IDH1),
    status = as.numeric(status)
  )

# Excluír patient_id se non se usa
data <- data %>% select(-patient_id)  # Opcional: Eliminar completamente

# --------------------------------------------
# PASO 3: División dos datos (20/80)
# --------------------------------------------
set.seed(123)
train_indices <- sample(1:nrow(data), 0.8 * nrow(data))
train_data <- data[train_indices, ]
test_data <- data[-train_indices, ]

# --------------------------------------------
# PASO 4: Entrenamento do modelo (CORRECCIÓN)
# --------------------------------------------
modelo <- rfsrc(
  Surv(survival_time, status) ~ .,  # Agora non inclúe patient_id
  data = train_data,
  ntree = 500,
  importance = TRUE,
  splitrule = "logrank"
)

# --------------------------------------------
# PASO 5: Análise de importancia de variables
# --------------------------------------------
# Importancia VIMP (Variable Importance)
vimp <- vimp(modelo, joint = FALSE)
print(vimp)

# Gráfico de importancia
ggplot(data = data.frame(Variable = names(vimp$importance), Importance = vimp$importance),
       aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Importancia de Variables no Modelo", x = "Variable", y = "Importancia (VIMP)") +
  theme_minimal()

ggsave("variable_importance.png", width = 10, height = 6)

# --------------------------------------------
# PASO 6: Predición e estratificación por risco
# --------------------------------------------
# Predición no conxunto de proba
predicion <- predict(modelo, newdata = test_data)

# Dividir en grupos de risco (alto/baixo)
risk_groups <- ifelse(predicion$predicted > median(predicion$predicted), "Alto Risco", "Baixo Risco")
test_data$risk_group <- factor(risk_groups)

# --------------------------------------------
# PASO 7: Curvas de supervivencia estratificadas (Kaplan-Meier)
# --------------------------------------------
# Axustar modelo Kaplan-Meier
km_fit <- survfit(Surv(survival_time, status) ~ risk_group, data = test_data)

# Graficar
ggsurvplot(
  km_fit,
  data = test_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  palette = c("#E7B800", "#2E9FDF"),
  title = "Curvas de Supervivencia por Grupo de Risco",
  xlab = "Tempo (meses)",
  ylab = "Probabilidade de Supervivencia"
)

ggsave("curvas_supervivencia.png", width = 10, height = 6)

# --------------------------------------------
# PASO 8: Validación do modelo (Opcional)
# --------------------------------------------
# C-index para medir desempeño
c_index <- rfsrc(Surv(survival_time, status) ~ ., data = test_data)$err.rate[500]
print(paste("C-index do modelo:", c_index))