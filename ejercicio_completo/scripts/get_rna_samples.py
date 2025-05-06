import os
import re

def get_rna_samples(data_dir):
    """
    Extrae IDs de muestra y tipos de tumor de archivos RNA FASTQ.
    Formato esperado: RXXX.RNA.<tumor_type>_N.fastq.gz
    """
    samples = set()
    
    for filename in os.listdir(data_dir):
        if not filename.endswith(".fastq.gz") or ".RNA." not in filename:
            continue
        
        # Extraer componentes usando expresión regular
        match = re.match(r"R(\d+)\.RNA\.(primary_tumor|recurrent_tumor)_\d\.fastq\.gz", filename)
        if match:
            sample_id = match.group(1)
            tumor_type = match.group(2)
            samples.add((sample_id, tumor_type))
    
    return sorted(samples, key=lambda x: (x[0], x[1]))

if __name__ == '__main__':
    data_dir = "/mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/exercicio_completo/data"
    samples = get_rna_samples(data_dir)
    print("Muestras de ARN encontradas:")
    for sample_id, tumor_type in samples:
        print(f"ID: {sample_id}, Tipo: {tumor_type}")