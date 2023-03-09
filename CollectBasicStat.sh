#!/bin/bash

# Script for collecting basic statistics of Mus Musculus mm9 genome

# Number of sequences in the reference genome:
NUM_SEQS=$(grep -c "^>" "/home/bioinformatikai/HW1/references/mm9.fa.gz")
echo "Number of sequences in the reference genome: ${NUM_SEQS}"

# Number of reads in each sample:
for FASTQ in /home/bioinformatikai/HW1/inputs/*.fastq.gz; do
  NUM_READS=$(zcat "${FASTQ}" | awk '{s++}END{print s/4}')
  echo "Number of reads in $(basename ${FASTQ}): ${NUM_READS}"
done 

# Number of protein-coding genes in the genome:
NUM_GENE=$(zcat /home/bioinformatikai/HW1/references/mm9.gtf.gz | awk '$3=="gene" && $0~/protein_coding/' | uniq -f 8 | wc -l)
echo "Number of of protein-coding genes in the genome: ${NUM_GENE}"
