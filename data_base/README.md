# Análisis de base de datos pública mediante R

En este directorio se encuentra una pequeña prueba de empleo de modelos de Machine Learning usando el lenguaje de programación R. 
Se emplea una base de datos de pacientes entre los cuales hay pacientes sanos y con diabetes. En esta base de datos se almacenan
una serie de parámetros que se obtienen de un análisis sanguíneo. El objetivo de este ejercicio es la obtención del (o los) parámetros
que diferencian el grupo control del grupo con enferemdad. 

Se ha separado el ejercicio en dos partes:

1. Análisis mediante RandomForest de la base de datos: se usa este modelo por ser uno de los más conocidos, a pesar de no ser el más
   adecuado (como se indica a continuación). Se obtiene así un Resultados_RandomForest_Diabetes_CON_Supervivencia.pdf que muestra los
   resultados obtenidos en este análisis, mostrando las diferencias entre un grupo y otro.
2. Comparación de modelos de Machine Learning: se comparan 3 modelos para el análisis de la base de datos: RandomForest, SVM y XGB.
   Se observa gran similitud entre los 3 modelos seleccionados, con una pequeña ventaja para XGB. Sin embargo, las diferencias no son
   significativas.
