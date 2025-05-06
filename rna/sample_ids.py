import os
from os.path import isfile, join

def get_sample_ids(data):
    unique_ids = set()
    paired_ids = set()

    for filename in os.listdir(data):
        if filename.endswith('_R1_001.fastq.gz') or filename.endswith('_R1.fastq.gz'):
            if '_R1_001.fastq.gz' in filename:
                file_id = filename.rsplit('_R1_001.fastq.gz', 1)[0]
            else:
                file_id = filename.rsplit('_R1.fastq.gz', 1)[0]
            unique_ids.add((file_id, 'R1'))
        elif filename.endswith('_R2_001.fastq.gz') or filename.endswith('_R2.fastq.gz'):
            if '_R2_001.fastq.gz' in filename:
                file_id = filename.rsplit('_R2_001.fastq.gz', 1)[0]
            else:
                file_id = filename.rsplit('_R2.fastq.gz', 1)[0]
            unique_ids.add((file_id, 'R2'))

    id_dict = {}
    for file_id, read_type in unique_ids:
        if file_id not in id_dict:
            id_dict[file_id] = set()
        id_dict[file_id].add(read_type)

    for file_id, read_types in id_dict.items():
        if 'R1' in read_types and 'R2' in read_types:
            paired_ids.add(file_id)

    sorted_paired_ids =sorted(list(paired_ids))
    return sorted_paired_ids        