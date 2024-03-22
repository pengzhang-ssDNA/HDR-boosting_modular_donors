import gzip
from Bio import SeqIO
import os

# Directory containing the input files
input_dir = 'D:/R/ISSD_Donor/NC/20240315_Tn5/FANCF'

# Directory for the output files
output_dir = 'D:/R/ISSD_Donor/NC/20240315_Tn5/FANCF/Filter'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Sequence to filter by (Read2 must start with this sequence)
sequence_to_filter = "CCCTTCTGCAGCAAGCTTATCGCTTTTCCGAG"

def read2_starts_with_sequence(record, sequence):
    """Check if Read2 starts with the sequence."""
    return str(record.seq).startswith(sequence)

def filter_reads(input_r1_path, input_r2_path, output_r1_path, output_r2_path):
    """Filter reads where Read2 starts with a specific sequence."""
    with gzip.open(input_r1_path, "rt") as input_r1, \
         gzip.open(input_r2_path, "rt") as input_r2, \
         gzip.open(output_r1_path, "wt") as output_r1, \
         gzip.open(output_r2_path, "wt") as output_r2:

        r1_records = SeqIO.parse(input_r1, "fastq")
        r2_records = SeqIO.parse(input_r2, "fastq")

        for r1_record, r2_record in zip(r1_records, r2_records):
            if read2_starts_with_sequence(r2_record, sequence_to_filter):
                SeqIO.write(r1_record, output_r1, "fastq")
                SeqIO.write(r2_record, output_r2, "fastq")

# Iterate over all file pairs in the input directory
for file in os.listdir(input_dir):
    if file.endswith("_1.fq.gz"):
        base_name = file[:-len("_1.fq.gz")]
        r1_path = os.path.join(input_dir, file)
        r2_path = os.path.join(input_dir, f"{base_name}_2.fq.gz")
        
        if os.path.exists(r2_path):  # Ensure the pair exists
            output_r1_path = os.path.join(output_dir, f"{base_name}_filtered_R1.fq.gz")
            output_r2_path = os.path.join(output_dir, f"{base_name}_filtered_R2.fq.gz")
            print(f"Filtering {file} and its pair...")
            filter_reads(r1_path, r2_path, output_r1_path, output_r2_path)

print("All filtering complete.")
