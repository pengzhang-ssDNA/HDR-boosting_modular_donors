#!/bin/bash

# Activate the conda environment with necessary tools
source activate atac

# Directories and paths
DATA_DIR="/home/issd/Rdata/issd_donor/nc/Donor_OT/data"
OUT_DIR="${DATA_DIR}/processed_bam"
BT2_INDEX="/home/issd/genomes/bowtie2_indexes/hg19"
REF_GENOME="/home/issd/genomes/hg19/hg19.fa"
mkdir -p "$OUT_DIR"

# Iterate over R1 files
for R1 in "${DATA_DIR}"/*_R1.fq.gz; do
    SAMPLE_NAME=$(basename "${R1}" "_R1.fq.gz")
    
    # Define output files
    SAM="${OUT_DIR}/${SAMPLE_NAME}.sam"
    BAM="${OUT_DIR}/${SAMPLE_NAME}.bam"
    SORTED_BAM="${OUT_DIR}/${SAMPLE_NAME}_sorted.bam"
    SORTED20_BAM="${OUT_DIR}/${SAMPLE_NAME}_sorted20.bam"
    REGION_BAM="${OUT_DIR}/${SAMPLE_NAME}_region.bam"
    DELLY_BCF="${OUT_DIR}/${SAMPLE_NAME}_translocations.bcf"
    DELLY_VCF="${OUT_DIR}/${SAMPLE_NAME}_translocations.vcf"
    REGION_TXT="${OUT_DIR}/${SAMPLE_NAME}_region.txt"
    source activate atac
    # Bowtie2: Map reads to the reference genome
    bowtie2 --very-sensitive -x "${BT2_INDEX}" -U "${R1}" -S "${SAM}"
    
    # SAMtools: Convert SAM to BAM, sort, and index
    samtools view -bS "${SAM}" > "${BAM}"
    samtools sort "${BAM}" -o "${SORTED_BAM}"
    samtools index "${SORTED_BAM}"
  
    # Filter unique reads with mapping length > 20
    samtools view -h "${SORTED_BAM}" | awk 'length($10) > 20 || $1 ~ /^@/' | samtools view -bS - > "${SORTED20_BAM}"
    
    samtools view -h "${SORTED20_BAM}" | awk '($1 !~ /^@/ && $6 !~ /H/) {chr=$3; start=$4; end=start+length($10); region=int(start/100000)*100000; print chr"\t"region"\t"region+100000"\t1"}' | sort -k1,1 -k2,2n | uniq -c | awk '{print $2"\t"$3"\t"$4"\t"$1}' > "${REGION_BAM}"
    
    samtools depth "${REGION_BAM}" > "${REGION_TXT}"
    # Optional: Remove intermediate SAM and BAM files
    rm "${SAM}" "${BAM}"
    
    # DELLY: Detect translocations
    delly call -t TRA -g "${REF_GENOME}" -o "${DELLY_BCF}" "${SORTED_BAM}"

    source activate ot
    
    # bcftools: Convert BCF to VCF for easier inspection
    bcftools view "${DELLY_BCF}" > "${DELLY_VCF}"
    
    echo "Translocations detected and processed for ${SAMPLE_NAME}"
done

echo "All samples processed."
