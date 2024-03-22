import csv
import gzip
import os
from Bio import SeqIO

folder_path = 'D:/R/ISSD_Donor/Insert/RUNX1'  # Replace with the actual folder path
target_sequence = 'GTGGGTACGAAGGAAATGACTCAA'

def count_reads(file_path, target_sequence):
    total_reads = 0
    double_occurrence_count = 0

    with gzip.open(file_path, 'rt') as fastq_file:
        for record in SeqIO.parse(fastq_file, 'fastq'):
            total_reads += 1
            seq = str(record.seq)
            if seq.count(target_sequence) == 2:
                double_occurrence_count += 1

    return total_reads, double_occurrence_count

# Process the folder
results = []
for filename in os.listdir(folder_path):
    if filename.endswith('.fq.gz'):  # Filter only compressed FASTQ files
        file_path = os.path.join(folder_path, filename)
        total_reads, double_occurrence_count = count_reads(file_path, target_sequence)
        proportion = double_occurrence_count / total_reads * 100 if total_reads != 0 else 0
        results.append([filename, total_reads, double_occurrence_count, f'{proportion:.2f}%'])

# Prepare the data to be written to the CSV file
data = [['File Name', 'Total Reads', 'Double Occurrence Count', 'Proportion']] + results

# Write the data to a CSV file
output_file = 'D:/R/ISSD_Donor/Insert/RUNX1/donor_insertion.csv'  # Replace with the desired output file path
with open(output_file, 'w', newline='') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerows(data)

print(f"Results saved successfully to: {output_file}")
