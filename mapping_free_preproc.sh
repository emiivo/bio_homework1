#!/bin/bash

REFERENCE_DIR="/home/bioinformatikai/HW1/references"
TRANSCRIPTOME_FILE="${REFERENCE_DIR}/mm9_RNA.fa.gz"
TRANSCRIPTOME_INDEX="${REFERENCE_DIR}/mm9_RNA_index"
INPUT_DIR="/home/bioinformatikai/HW1/inputs"
OUTPUT_DIR="/home/bioinformatikai/HW1/outputs"
FILES=(SRR8985047 SRR8985048 SRR8985051 SRR8985052)

# CHeck if index exists for transcriptome, and if it does not, create one.

if [ ! -d "${TRANSCRIPTOME_INDEX}" ]
then
  salmon index --index "${TRANSCRIPTOME_INDEX}" --transcripts "${TRANSCRIPTOME_FILE}"
fi

# Run FASTQC analysis on each of your FASTQC files.
fastqc /home/bioinformatikai/HW1/inputs/*.fastq.gz

# Generate MULTIQC report for FASTQC analysis results.
multiqc /home/bioinformatikai/HW1/inputs/ -o /home/bioinformatikai/HW1/inputs/

# Run standard FASTQ trimming: remove adapters, trim low-quality bases as well as remove reads that are shorter than 20 bp. 
for file in "${FILES[@]}"
do
  # Define input and output file names
  INPUT1="${INPUT_DIR}/${file}_1.fastq.gz"
  INPUT2="${INPUT_DIR}/${file}_2.fastq.gz"
  OUTPUT1="${OUTPUT_DIR}/${file}_1_trimmed.fastq.gz"
  OUTPUT2="${OUTPUT_DIR}/${file}_2_trimmed.fastq.gz"

  # Run trim_galore
  trim_galore --length 20 --paired ${INPUT1} ${INPUT2}

  # Rename output files
  #mv ${file}_1_val_1.fq.gz ${OUTPUT1}
  #mv ${file}_2_val_2.fq.gz ${OUTPUT2}

  # Quantify samples
  salmon quant -i ${TRANSCRIPTOME_INDEX} -l A -1 ${OUTPUT1} -2 ${OUTPUT2} -o ${OUTPUT}
done

# Rerun FASTQC on newly created/cleaned FASTQ files.
fastqc /home/bioinformatikai/HW1/outputs/*.fastq.gz

#Create MultiQC plot for processed data
multiqc /home/bioinformatikai/HW1/outputs/ -o /home/bioinformatikai/HW1/outputs/
