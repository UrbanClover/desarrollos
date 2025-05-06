import pandas as pd
import numpy as np

# Configuración inicial
np.random.seed(42)
n_pacientes = 5000

# 1. ID do paciente
ids = [f"P{str(i).zfill(4)}" for i in range(1, n_pacientes + 1)]

# 2. Idade e sexo
idade = np.clip(np.random.normal(60, 15, n_pacientes).astype(int), 18, 90)
sexo = np.random.choice([0, 1], n_pacientes, p=[0.5, 0.5])  # 0: Muller, 1: Home

# 3. Mutación en LTBP4/IDH1 (30% con mutación)
mutacion = np.random.binomial(1, 0.3, n_pacientes)

# 4. Nivel de expresión de LTBP4 (log-normal)
ltbp4_exp = np.round(np.random.lognormal(mean=3, sigma=1, size=n_pacientes))

# 5. Expresión dos 10 transcritos (modelado baseado en datos orixinais)
transcricions = [
    "ENST00000410691.1", "ENST00000416931.1", "ENST00000457540.1",
    "ENST00000514057.1", "ENST00000523172.5", "ENST00000579103.1",
    "ENST00000619216.1", "ENST00000620525.1", "ENST00000621981.1",
    "ENST00000623083.4"
]

# Exemplo de modelo para un transcrito (ENST00000410691.1 con 50% ceros)
def simular_transcrito(prob_cero, media_non_cero, var_non_cero, n):
    ceros = np.random.binomial(1, prob_cero, n)
    non_cero = np.round(np.random.gamma((media_non_cero**2)/var_non_cero, var_non_cero/media_non_cero, n))
    return np.where(ceros, 0, non_cero)

expresion_transcritos = {
    transcricions[0]: simular_transcrito(0.5, 6.3, 22, n_pacientes),
    # Engadir aquí os parámetros para os outros 9 transcritos
}

# 6. Variables clínicas simuladas
ldh = np.random.gamma(2, 150, n_pacientes).astype(int)  # U/L
glucosa = np.clip(np.random.normal(90, 15, n_pacientes), 70, 200).astype(int)  # mg/dL
hemoglobina = np.where(sexo == 1, 
                      np.random.normal(15, 1, n_pacientes),
                      np.random.normal(13, 1, n_pacientes)).round(1)  # g/dL
plaquetas = np.clip(np.random.normal(250, 50, n_pacientes), 150, 450).astype(int)  # 10^3/µL

# 7. Status e SurvivalTime
status = np.random.binomial(1, 0.3, n_pacientes)
tempo_supervivencia = np.where(
    status == 1,
    np.random.weibull(1.5, n_pacientes) * 30,  # Eventos
    np.random.uniform(50, 120, n_pacientes)     # Censura
).astype(int)

# Crear DataFrame
df = pd.DataFrame({
    "patient_id": ids,
    "idade": idade,
    "sexo": sexo,
    "mutacion_LTBP4_IDH1": mutacion,
    "ltbp4_expresion": ltbp4_exp,
    **expresion_transcritos,
    "LDH": ldh,
    "Glucosa": glucosa,
    "Hemoglobina": hemoglobina,
    "Plaquetas": plaquetas,
    "status": status,
    "survival_time": tempo_supervivencia
})

# Gardar en TSV
df.to_csv("simulated_patients.tsv", sep="\t", index=False)