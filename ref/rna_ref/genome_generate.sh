#!/bin/bash
#SBATCH --job-name=genome_index_primary_assembly
#SBATCH --output=genome_index_primary_assembly_%j.out
#SBATCH --error=genome_index_primary_assembly_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# Configuración de rutas y parámetros
STAR_INDEX_DIR="ref_index"
GENOME_FASTA="GRCh38.primary_assembly.genome.fa"
ANNOTATION_GTF="gencode.v47.primary_assembly.basic.annotation.gtf"
THREADS=256                         # Número de hilos para procesamiento paralelo
SJDB_OVERHANG=100                 # Longitud de lectura - 1 (ajustar según experimento)
GENOME_SA_SPARSE_D=2              # Para genomas grandes, valor 1 para pequeños
GENOME_SA_INDEX_NBASES=14
GENOME_CHR_BIN_NBITS=18

# Ejecutar comando STAR
module load cesga/2022 gcc/system star/2.7.11b
STAR --runThreadN $THREADS \
     --runMode genomeGenerate \
     --genomeDir "$STAR_INDEX_DIR" \
     --genomeFastaFiles "$GENOME_FASTA" \
     --sjdbGTFfile "$ANNOTATION_GTF" \
     --sjdbOverhang $SJDB_OVERHANG \
     --genomeSAsparseD $GENOME_SA_SPARSE_D \
     --genomeSAindexNbases $GENOME_SA_INDEX_NBASES \
     --genomeChrBinNbits $GENOME_CHR_BIN_NBITS