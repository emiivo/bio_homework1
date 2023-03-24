#!/usr/bin

# Check if index files exist and if it does not exist - create a reference genome index. 
REFERENCE_DIR="/home/bioinformatikai/HW1/references"
GENOME_FILE="${REFERENCE_DIR}/mm9.fa.gz"
GENOME_UNZIPPED="${REFERENCE_DIR}/mm9.fa"
GENOME_INDEX="${REFERENCE_DIR}/mm9"
INPUT_DIR="/home/bioinformatikai/HW1/inputs"
OUTPUT_DIR="/home/bioinformatikai/HW1/outputs"
GTF_FILE="/home/bioinformatikai/HW1/references/mm9.gtf.gz"
GTF_FILE_UNZIPPED="/home/bioinformatikai/HW1/references/mm9.gtf"
FILES=(SRR8985047 SRR8985048 SRR8985051 SRR8985052)

gunzip -c "${GTF_FILE}" > "${GTF_FILE_UNZIPPED}"

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
  OUTPUT_SAM_FILE="${OUTPUT_DIR}/${file}.sam"
  INPUT_FILE_1="${OUTPUT_DIR}/${file}_1_trimmed.fastq.gz"
  INPUT_FILE_2="${OUTPUT_DIR}/${file}_2_trimmed.fastq.gz"
  INPUT_SAM_FILE="${OUTPUT_DIR}/${file}.sam"
  OUTPUT_BAM_FILE="${OUTPUT_DIR}/${file}.bam"

  # Run HISAT2
  hisat2 -p 6 --dta -x "${GENOME_INDEX}" -1 "$INPUT_FILE_1" -2 "$INPUT_FILE_2" -S "$OUTPUT_SAM_FILE"
  
  # Convert SAM to BAM format, remove unmapped reads as well as reads that are mapped incorrectly
  samtools view -@ 6 -F 0x4 -F 0x2 -bS "$INPUT_SAM_FILE" > "$OUTPUT_BAM_FILE"
done

# Deduplicate and index BAM files
for file in "${FILES[@]}"
do
  # Define input and output file names
  INPUT_BAM_FILE="${OUTPUT_DIR}/${file}.bam"
  SORTED_BAM_FILE="${OUTPUT_DIR}/${file}_sorted.bam"
  DEDUPLICATED_BAM_FILE="${OUTPUT_DIR}/${file}_deduplicated.bam"
  INDEX_FILE="${OUTPUT_DIR}/${file}_deduplicated.bam.bai"
  FIXED_BAM_FILE="${OUTPUT_DIR}/${file}_fixed.bam"
  SORTED_BAM_FILE_2="${OUTPUT_DIR}/${file}_sorted_2.bam"
 
  
  samtools sort -n "${INPUT_BAM_FILE}" -o "${SORTED_BAM_FILE}"
  samtools fixmate -m "${SORTED_BAM_FILE}" "${FIXED_BAM_FILE}"
  samtools sort "${FIXED_BAM_FILE}" -o "${SORTED_BAM_FILE_2}"
  samtools markdup -r "${SORTED_BAM_FILE_2}" "${DEDUPLICATED_BAM_FILE}"
  # Index the deduplicated BAM file
  samtools index -b "${DEDUPLICATED_BAM_FILE}" "${INDEX_FILE}"
  
  rm "${SORTED_BAM_FILE}" "${FIXED_BAM_FILE}"
  
  # Run StringTie
  stringtie -G "${GTF_FILE_UNZIPPED}" -o "${OUTPUT_DIR}/${file}.gtf" "${SORTED_BAM_FILE_2}"

done
# Run multiBamSummary
multiBamSummary bins --outFileName "${RESULTS_DIR}/mapped.npz" --binSize 1000 -p 6 --outRawCounts "${RESULTS_DIR}/raw_counts.tsv" -b "${OUTPUT_DIR}/SRR8985047_deduplicated.bam" "${OUTPUT_DIR}/SRR8985048_deduplicated.bam" "${OUTPUT_DIR}/SRR8985051_deduplicated.bam" "${OUTPUT_DIR}/SRR8985052_deduplicated.bam"

# Generate the correlation plots
plotCorrelation -in "${RESULTS_DIR}/mapped.npz"  --plotNumbers -c pearson -p heatmap -o "${RESULTS_DIR}/mapped_data_heatmap.pdf"
plotCorrelation -in "${RESULTS_DIR}/mapped.npz" -c pearson -p scatterplot -o "${RESULTS_DIR}/mapped_data_scatterplot.pdf"

# Generate the PCA plot
plotPCA -in "${RESULTS_DIR}/mapped.npz" -o "${RESULTS_DIR}/mapped_data_PCA.pdf"
