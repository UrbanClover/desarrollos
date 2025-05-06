import os

def get_sample_ids(data):
    """
    Extrae los IDs de muestra de los archivos de ADN (DNA) en el directorio 'data'.
    Solo se consideran aquellos archivos que contengan '.DNA.' y que tengan finalización
    '_1.fastq.gz' o '_2.fastq.gz'. Se emparejan los archivos de lectura 1 (R1) y lectura 2 (R2)
    y se retorna una lista ordenada con los IDs que cuenten con ambos pares.
    """
    unique_ids = set()
    paired_ids = set()

    for filename in os.listdir(data):
        # Procesar solo archivos de ADN
        if ".DNA." not in filename:
            continue

        # Verificar las lecturas R1 y R2 basadas en el sufijo de la nomenclatura
        if filename.endswith("_1.fastq.gz"):
            file_id = filename.rsplit("_1.fastq.gz", 1)[0]
            unique_ids.add((file_id, "R1"))
        elif filename.endswith("_2.fastq.gz"):
            file_id = filename.rsplit("_2.fastq.gz", 1)[0]
            unique_ids.add((file_id, "R2"))

    # Organizar los IDs para verificar la existencia del par de archivos (R1 y R2)
    id_dict = {}
    for file_id, read_type in unique_ids:
        if file_id not in id_dict:
            id_dict[file_id] = set()
        id_dict[file_id].add(read_type)

    # Solo se consideran aquellos IDs que tengan ambos tipos de lectura
    for file_id, read_types in id_dict.items():
        if "R1" in read_types and "R2" in read_types:
            paired_ids.add(file_id)

    # Retornar la lista de IDs emparejados ordenada
    sorted_paired_ids = sorted(list(paired_ids))
    return sorted_paired_ids

if __name__ == '__main__':
    # Ejemplo de uso: se asume que 'data' es el directorio donde se encuentran los archivos .fastq.gz
    data_dir = "/mnt/lustre/scratch/nlsas/home/otras/hcx/iba/ngs/exercicio_completo/data"
    samples = get_sample_ids(data_dir)
    print("Muestras de ADN encontradas:")
    for sample in samples:
        print(sample)
