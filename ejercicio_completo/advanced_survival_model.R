# --------------------------------------------
# PASO 1: Instalar e cargar paquetes (CORRECCIÓN FINAL)
# --------------------------------------------
if (!require("randomForestSRC")) install.packages("randomForestSRC")
if (!require("survival")) install.packages("survival")
if (!require("survminer")) install.packages("survminer")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggpubr")) install.packages("ggpubr")
if (!require("cowplot")) install.packages("cowplot")  # ESSENCIAL PARA ggdraw()
if (!require("magick")) install.packages("magick")

library(randomForestSRC)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot) 
library(magick)  

# --------------------------------------------
# PASO 2: Carga e Preparación de Datos
# --------------------------------------------
data <- read.delim("simulated_patients.tsv", sep = "\t", stringsAsFactors = TRUE) %>%
  mutate(
    sexo = as.factor(sexo),
    mutacion_LTBP4_IDH1 = as.factor(mutacion_LTBP4_IDH1),
    status = as.numeric(status)
  ) %>%
  select(-patient_id) %>%
  na.omit()

# --------------------------------------------
# PASO 3: División Estratificada (80/20)
# --------------------------------------------
set.seed(123)
train_indices <- sample(1:nrow(data), 0.8 * nrow(data))
train_data <- data[train_indices, ]
test_data <- data[-train_indices, ]

# --------------------------------------------
# PASO 4: Tuning de Hiperparámetros
# --------------------------------------------
tune_results <- tune(
  Surv(survival_time, status) ~ .,
  data = train_data,
  ntreeTry = 500,
  mtryStart = 5,
  nodesizeTry = c(5, 10, 20),
  doBest = FALSE
)

# --------------------------------------------
# PASO 5: Entrenamento do Modelo Final
# --------------------------------------------
modelo_final <- rfsrc(
  Surv(survival_time, status) ~ .,
  data = train_data,
  ntree = 1000,
  mtry = tune_results$optimal[["mtry"]],
  nodesize = tune_results$optimal[["nodesize"]],
  splitrule = "logrank",
  importance = TRUE
)

# --------------------------------------------
# PASO 6: Análise de Importancia
# --------------------------------------------
vimp <- vimp(modelo_final, joint = FALSE)

ggplot(data.frame(Variable = names(vimp$importance), Importance = vimp$importance),
       aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Importancia de Variables", x = "Variable", y = "VIMP") +
  theme_minimal()

ggsave("variable_importance.png", width = 10, height = 6)

# --------------------------------------------
# PASO 7: Predición e Curvas de Supervivencia
# --------------------------------------------
predicion <- predict(modelo_final, newdata = test_data)
test_data$predicted_risk <- predicion$predicted

# Estratificación en 2 grupos
test_data$risk_group <- ifelse(test_data$predicted_risk > median(test_data$predicted_risk), 
                              "Alto Risco", "Baixo Risco")

# Kaplan-Meier
km_fit <- survfit(Surv(survival_time, status) ~ risk_group, data = test_data)

ggsurvplot(
  km_fit,
  data = test_data,
  pval = TRUE,
  conf.int = TRUE,
  palette = c("#E7B800", "#2E9FDF"),
  title = "Curvas de Supervivencia por Risco",
  xlab = "Meses"
)

ggsave("curvas_supervivencia.png", width = 10, height = 6)

# --------------------------------------------
# PASO 8: Gráfico de Converxencia do Error
# --------------------------------------------
error_df <- data.frame(
  Arbores = 1:modelo_final$ntree,
  Error = modelo_final$err.rate
)

ggplot(error_df, aes(x = Arbores, y = Error)) +
  geom_line(color = "steelblue") +
  labs(title = "Converxencia do Error do Modelo", x = "Número de Árbores", y = "Error") +
  theme_minimal()

ggsave("error_convergence.png", width = 10, height = 6)

# --------------------------------------------
# PASO 9: Validación con C-index
# --------------------------------------------
c_index <- 1 - modelo_final$err.rate[length(modelo_final$err.rate)]
print(paste("C-index:", round(c_index, 3)))

# --------------------------------------------
# PASO 10: Informe Final en PDF
# --------------------------------------------
pdf("Informe_Final_Modelo.pdf", width = 12, height = 8)

# Páxina 1: Importancia de Variables
p1 <- ggdraw() + 
  draw_image("variable_importance.png") +
  draw_label("Análise de Importancia de Variables", x = 0.5, y = 0.95, size = 14)

# Páxina 2: Curvas de Supervivencia + C-index
p2 <- ggdraw() + 
  draw_image("curvas_supervivencia.png") +
  draw_label(paste("C-index:", round(c_index, 3)), x = 0.5, y = 0.05, size = 12)

# Páxina 3: Converxencia do Error
p3 <- ggdraw() + 
  draw_image("error_convergence.png")

print(p1)
print(p2)
print(p3)

dev.off()