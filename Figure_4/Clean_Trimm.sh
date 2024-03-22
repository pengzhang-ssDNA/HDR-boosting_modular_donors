#!/bin/bash
# Set the path to Trimmomatic
trimmomatic_path="/home/issd/anaconda3/envs/crispr/share/trimmomatic-0.39-2/trimmomatic"
# Set the input and output folder paths
input_folder="$HOME/Rdata/ISSD_Donor/nc/20240315_Tn5/Raw_data"
output_folder="$HOME/Rdata/ISSD_Donor/nc/20240315_Tn5/Clean_data"
# Create the output folder if it doesn't exist
mkdir -p "$output_folder"
# Function to handle file processing for paired-end data
process_files() {
  local file_r1="$1"
  local file_r2="$2"
  local base_name=$(basename -- "$file_r1" | cut -d "_" -f 1-5)
  base_name="${base_name%_1.fq.gz}"

  local output_file_r1="$output_folder/${base_name}_1.fq.gz"
  local output_file_r2="$output_folder/${base_name}_2.fq.gz"
  local output_file_r1_unpaired="$output_folder/${base_name}_R1_clean_unpaired.fq.gz"
  local output_file_r2_unpaired="$output_folder/${base_name}_R2_clean_unpaired.fq.gz"

  # Run Trimmomatic for paired-end data
  trimmomatic PE -threads 12 \
    "$file_r1" "$file_r2" \
    "$output_file_r1" "$output_file_r1_unpaired" \
    "$output_file_r2" "$output_file_r2_unpaired" \
    ILLUMINACLIP:adapter.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

  if [[ $? -eq 0 ]]; then
    echo "Processed: $file_r1 and $file_r2 -> $output_file_r1, $output_file_r2"
  else
    echo "Error: Trimmomatic failed on '$file_r1' and '$file_r2'."
    return 1
  fi
}

# Loop through all pairs of .fq.gz files in the input folder
for file_r1 in "$input_folder"/*_1.fq.gz; do
  file_r2="${file_r1/_1/_2}"
  if [[ -f "$file_r2" ]]; then
    process_files "$file_r1" "$file_r2"
  else
    echo "Error: Paired file for '$file_r1' not found."
  fi
done

echo "Trimming complete!"

exit 0
