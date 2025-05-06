# Desde el directorio del proyecto:
Rscript scripts/differential_analysis.R \
  -i /mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/exercicio_completo/results/expression_matrix/transcript_matrix.tsv \
  -m /mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/exercicio_completo/results/expression_matrix/sample_metadata.tsv \
  -o results/top_transcripts/top10_common_transcripts.txt \
  -d "~ tumor_type"