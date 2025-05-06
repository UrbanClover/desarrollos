# Análisis de Muestras de ARN

En este directorio se encuentran los ficheros para el análisis de muestras de ARN.

En data, concretamente en descarga_muestras.sh, están enlazados para su descarga las muestras que se han usado en este ejemplo.

En el archivo Snakefile se encuentra el código para el procesado de las muestras de este ejemplo. Para una correcta indentificación
de las mismas, es necesario que el archivo sample_ids.py se encuentre en el directorio rna. Con el archivo Snakefile se realizan los
trabajos relacionados con STAR, Salmon y Arriba. 

Posteriormente, con los archivos resultantes de estos análisis, y mediante los diversos scripts de R que se encuentran en este directorio,
se generan una serie de archivos y gráficos (encuadrados tanto en la carpeta principal como en la carpeta Scripts). Estos archivos y 
gráficos finales muestran los resultados de comparar las muestras control con las muestras de cáncer, mostrando los puntos claves que
diferencian ambos grupos y los genes implicados.
