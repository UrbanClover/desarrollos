#!/bin/bash
#SBATCH --job-name=indices_rna
#SBATCH --output=indices_rna_%j.out
#SBATCH --error=indices_rna_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# -----------------------------
# Definición de variables
# -----------------------------
# Archivo FASTA con el transcriptoma (para Salmon)
TRANSCRIPTOMA="GRCh38.p14.genome.fa"

# Directorio donde se creará el índice de Salmon
SALMON_INDEX="salmon_index/"

# Archivo FASTA del genoma (para Arriba)
GENOMA="GRCh38.p14.genome.fa"

# Archivo GTF de anotación (para Arriba)
ANOTACION="gencode.v47.chr_patch_hapl_scaff.annotation.gtf"

# Archivo de salida para el índice de Arriba
ARRIBA_INDEX="arriba/arriba_index.idx"

# -----------------------------
# Creación del índice para Salmon
# -----------------------------
echo "Creando índice para Salmon..."

# Si no usas decoys:
module load cesga/2020 gcc/system salmon/1.10.2
salmon index -t "$TRANSCRIPTOMA" -i "$SALMON_INDEX" -k 31
module unload cesga/2020 gcc/system salmon/1.10.2

# Si deseas usar decoys (descomenta la línea siguiente y comenta la anterior)
# salmon index -t "$TRANSCRIPTOMA" -i "$SALMON_INDEX" -d "$DECOYS" --gencode

# -----------------------------
# Creación del índice para Arriba
# -----------------------------
echo "Creando índice para Arriba..."

# Se asume que Arriba incluye un script de indexación llamado 'arriba_index.py'
# Asegúrate de que este script esté en el PATH o especifica su ruta completa.
module load cesga/2020 gcc/system arriba/2.4.0
python3 /ruta/al/script/arriba_index.py -g "$GENOMA" -a "$ANOTACION" -o "$ARRIBA_INDEX"
module unload cesga/2020 gcc/system arriba/2.4.0

echo "Proceso completado: índices creados correctamente."
