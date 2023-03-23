#!/usr/bin/env bash

# Check if index files exist and if it does not exist - create a reference genome index. 
REFERENCE_DIR="/home/bioinformatikai/HW1/references"
GENOME_FILE="${REFERENCE_DIR}/mm9.fa.gz"
GENOME_UNZIPPED="${REFERENCE_DIR}/mm9.fa"
GENOME_INDEX="${REFERENCE_DIR}/mm9"


if [ ! -f "${GENOME_UNZIPPED}" ]; then
    gunzip -c ${GENOME_FILE} > ${GENOME_UNZIPPED}
fi

# Create HISAT2 index if it doesn't exist
if [ ! -f "${GENOME_INDEX}.1.ht2" ]; then
    hisat2-build ${GENOME_UNZIPPED} ${GENOME_INDEX}
fi


# Run FASTQC analysis on each of your FASTQC files.
fastqc /home/bioinformatikai/HW1/inputs/*.fastq.gz

# Generate MULTIQC report for FASTQC analysis results.
multiqc /home/bioinformatikai/HW1/inputs/ -o /home/bioinformatikai/HW1/inputs/

# Run standard FASTQ trimming: remove adapters, trim low-quality bases as well as remove reads that are shorter than 20 bp. 
# Define input directory
INPUT_DIR=/home/bioinformatikai/HW1/inputs
OUTPUT_DIR=/home/bioinformatikai/HW1/outputs

# Define input file names
FILES=(SRR8985047 SRR8985048 SRR8985051 SRR8985052)

# Loop through each file and perform trimming
for file in "${FILES[@]}"
do
  # Define input and output file names


  # Run trim_galore
  trim_galore --length 20 --paired ${INPUT1} ${INPUT2}

  # Rename output files
  mv ${file}_1_val_1.fq.gz ${OUTPUT1}
  mv ${file}_2_val_2.fq.gz ${OUTPUT2}
done

# Rerun FASTQC on newly created/cleaned FASTQ files.
fastqc /home/bioinformatikai/HW1/outputs/*.fastq.gz

#Create MultiQC plot for processed data
multiqc /home/bioinformatikai/HW1/outputs/ -o /home/bioinformatikai/HW1/outputs/

# Mapping using HISAT2
for file in "${FILES[@]}"
do
  # Define input and output file names for hisat2
  OUTPUT_SAM_FILE="${OUTPUT_DIR}/${file}.sam"
  INPUT_FILE_1="${OUTPUT_DIR}/${file}_1_trimmed.fastq.gz"
  INPUT_FILE_2="${OUTPUT_DIR}/${file}_2_trimmed.fastq.gz"

  # Run HISAT2
  hisat2 -p 6 --dta -x "${GENOME_INDEX}" -1 "$INPUT_FILE_1" -2 "$INPUT_FILE_2" -S "$OUTPUT_SAM_FILE"
  
  # Define input and output file names for samtools
  INPUT_SAM_FILE="${OUTPUT_DIR}/${file}.sam"
  OUTPUT_BAM_FILE="${OUTPUT_DIR}/${file}.bam"

  # Convert SAM to BAM format
  samtools view -@ 6 -F 0x4 -F 0x2 -bS "$INPUT_SAM_FILE" > "$OUTPUT_BAM_FILE"
done
