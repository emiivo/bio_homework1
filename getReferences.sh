#!/usr/bin
 
#Download data of Mus Musculus version mm9

#Download the reference genome in FASTA format:
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.fa.gz -P /home/bioinformatikai/HW1/references/

#Download the reference transcriptome in FASTA format:
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/refMrna.fa.gz -O /home/bioinformatikai/HW1/references/mm9_RNA.fa.gz

#Download GTF/GFF3 file for your reference genome:
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/genes/mm9.refGene.gtf.gz -P /home/bioinformatikai/HW1/references/

#Download raw FASTQ files:
prefetch -O /home/bioinformatikai/HW1/inputs SRR8985047 SRR8985048 SRR8985051 SRR8985052

fastq-dump --outdir /home/bioinformatikai/HW1/inputs/ --gzip --split-files /home/bioinformatikai/HW1/inputs/SRR8985047 /home/bioinformatikai/HW1/inputs/SRR8985048 /home/bioinformatikai/HW1/inputs/SRR8985051 /home/bioinformatikai/HW1/inputs/SRR8985052


