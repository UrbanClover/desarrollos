# Análisis de Muestras de ADN

En este directorio se encuentran los ficheros para el análisis de muestras de ADN. 

En data, concretamente en descarga_muestras.sh, están enlazados para su descarga las muestras que se han usado en este ejemplo. 

Cada archivo de Snakemake, ordenado numéricamente por ejecución, tiene un archivo de bash que se ha empleado para su ejecución mediante mediante slurm en un servidor externo (CESGA en este caso).
Dichos archivos pueden usarse a modo de ejemplo para los parámetros de ejecución (pudiendo adaptarlos a cada caso). Se ha implementado esta división del trabajo en diversos archivos debido a que
hay ciertos programos (caso de Mutect2), que dan problemas en ejecuciones sucesivas, rompiendo el flujo de trabajo por errores. De este modo, se soluciona este problema de manera sencilla. También
se divide el trabajo para poder identificar problemas en el código de manera más sencilla.

En este caso, al ser una comparativa de muestras de ADN control y con cáncer, se realiza un procesado de las muestras y, posteriormente, se generan una serie de gráficos y archivos de salida que se
encuentran almacenados en los directorios resultados_annovar_multianno y scripts. Si se ejecutan los archivos de Snakemake que aquí se muestra, se generaran otros directorios de salida intermedios,
necesarios para la obtención de los archivos de finales que se han subido a este repositorio. 
