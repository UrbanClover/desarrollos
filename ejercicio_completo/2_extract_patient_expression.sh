cut -d "|" -f1 /mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/exercicio_completo/results/top_transcripts/top10_common_transcripts.txt > /mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/exercicio_completo/results/top_transcripts/top10_clean.txt

Rscript scripts/extract_patient_expression.R \
  -m /mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/exercicio_completo/results/expression_matrix/transcript_matrix.tsv \
  -t /mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/exercicio_completo/results/top_transcripts/top10_clean.txt \
  -s /mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/exercicio_completo/results/expression_matrix/sample_metadata.tsv \
  -o /mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/exercicio_completo/results/patient_expression/patient_expression_values.tsv