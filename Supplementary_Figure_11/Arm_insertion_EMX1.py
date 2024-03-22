import pandas as pd
from Bio import SeqIO
import os
import gzip
import re

folder_path = 'D:/R/ISSD_Donor/Insert/EMX1'
sequence = 'CCCCATTG.{1,24}GCCTGCTTCG' 
# "CCCCATTG" is the genomic sequence adjacent to the 5` homology arm fo donor; 
# "GCCTGCTTCG" is the terminal sequence of 5` homology arm that is incorporated with the modular sequence.

def process_fastq_file(file_path):
    total_reads = 0
    reads_with_sequence = 0
    matching_sequences = []  # List to store matching sequences

    with gzip.open(file_path, 'rt') as fastq_file:
        for record in SeqIO.parse(fastq_file, 'fastq'):
            total_reads += 1
            if re.search(sequence, str(record.seq)):
                reads_with_sequence += 1
                matching_sequences.append(str(record.seq))

    return total_reads, reads_with_sequence, matching_sequences

def process_folder(folder_path):
    results = []

    for filename in os.listdir(folder_path):
        if filename.endswith('.fq.gz'):
            file_path = os.path.join(folder_path, filename)
            result = process_fastq_file(file_path)
            results.append((filename, *result))

            # Save matching sequences to a separate CSV file
            matching_sequences_df = pd.DataFrame(result[2], columns=['Matching Sequences'])
            matching_sequences_file = os.path.join(folder_path, f'{filename}_matching_sequences.csv')
            matching_sequences_df.to_csv(matching_sequences_file, index=False)

    return results


results = process_folder(folder_path)
results_df = pd.DataFrame(results, columns=['File Name', 'Total Reads', f"Reads with sequence '{sequence}'", 'Percentage of reads with sequence'])

results_df['Percentage of reads with sequence'] = (results_df[f"Reads with sequence '{sequence}'"] / results_df['Total Reads']) * 100

output_file = 'D:/R/ISSD_Donor/Insert/EMX1/EMX1_arm_insert.csv'
results_df.to_csv(output_file, index=False)

print("Results saved successfully as CSV files.")