#!/bin/bash
OUTPUT_DIR=~/HW1/outputs
RESULTS_DIR=~/HW1/results
GTF_FILE_UNZIPPED=~/HW1/references/mm9.gtf
FILES=(SRR8985047 SRR8985048 SRR8985051 SRR8985052)



for file in "${FILES[@]}"
do

  SORTED_BAM_FILE_2="${OUTPUT_DIR}/${file}_sorted_2.bam"

  # Run StringTie with -e option to later generate counts with prepDE
  stringtie -G "${GTF_FILE_UNZIPPED}" -e -o "${OUTPUT_DIR}/${file}_forPrepDE_.gtf" "${SORTED_BAM_FILE_2}"

done
mkdir ~/HW1/outputs/gtf_files_prepDE

mv ~/HW1/outputs/*_forPrepDE.gtf ~/HW1/outputs/gtf_files_prepDE/

cd ~/HW1/outputs

python3 prepDE.py3 *_forPrepDE.gtf

sed -i 's/gtf_files_prepDE/SRR8985047/g; s/gtf_files_prepDE/SRR8985048/g; s/gtf_files_prepDE/SRR8985051/g; s/gtf_files_prepDE/SRR8985052/g' gene_count_matrix.csv
#!/bin/bash
OUTPUT_DIR=~/HW1/outputs
RESULTS_DIR=~/HW1/results
GTF_FILE_UNZIPPED=~/HW1/references/mm9.gtf
FILES=(SRR8985047 SRR8985048 SRR8985051 SRR8985052)



for file in "${FILES[@]}"
do

  SORTED_BAM_FILE_2="${OUTPUT_DIR}/${file}_sorted_2.bam"

  # Run StringTie with -e option to later generate counts with prepDE
  stringtie -G "${GTF_FILE_UNZIPPED}" -e -o "${OUTPUT_DIR}/${file}_forPrepDE_.gtf" "${SORTED_BAM_FILE_2}"

done
mkdir ~/HW1/outputs/gtf_files_prepDE

mv ~/HW1/outputs/*_forPrepDE.gtf ~/HW1/outputs/gtf_files_prepDE/

cd ~/HW1/outputs

python3 prepDE.py3 *_forPrepDE.gtf

sed -i 's/gtf_files_prepDE/SRR8985047/g; s/gtf_files_prepDE/SRR8985048/g; s/gtf_files_prepDE/SRR8985051/g; s/gtf_files_prepDE/SRR8985052/g' gene_count_matrix.csv
mv 

python3 prepDE.py3 [options] input_file.gtf

# Generate count matrix using tximport with salmon
Rscript -e '
  library(tximport)
  files <- c("'${OUTPUT_DIR}/SRR8985047_quantifyed/quant.sf'"
            ,"'${OUTPUT_DIR}/SRR8985048_quantifyed/quant.sf'"
            ,"'${OUTPUT_DIR}/SRR8985051_quantifyed/quant.sf'"
            ,"'${OUTPUT_DIR}/SRR8985052_quantifyed/quant.sf'")
            
  txi <- tximport(files, type="salmon", txOut=TRUE, ignoreTxVersion=TRUE, countsFromAbundance="lengthScaledTPM")
  write.table(txi$counts, file="'${RESULTS_DIR}/non_mapping_counts.tsv'"
             ,quote=FALSE,sep="\t",row.names=TRUE,col.names=NA)'

#fix count files
awk -F'\t' 'BEGIN {OFS=FS} NR>1 {split($1,a,"|"); $1=a[2]} 1' ${RESULTS_DIR}/non_mapping_counts.tsv > ${RESULTS_DIR}/correct_non_mapping_counts.tsv
sed -i '1s/.*/gene_id\tSRR8985047\tSRR8985048\tSRR8985051\tSRR8985052/' ${RESULTS_DIR}/correct_non_mapping_counts.tsv


awk -F',' 'BEGIN {OFS=FS} NR>1 {split($1,a,"|"); $1=a[1]} 1' ${OUTPUT_DIR}/gene_count_matrix.csv > ${OUTPUT_DIR}/correct_gene_count_matrix.csv
sed -i '1s/.*/gene_id,SRR8985047,SRR8985048,SRR8985051,SRR8985052/' ${OUTPUT_DIR}/correct_gene_count_matrix.csv
