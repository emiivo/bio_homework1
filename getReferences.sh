#!/usr/bin
 
#Download data of Mus Musculus version mm9

#Download the reference genome in FASTA format:
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M1/NCBIM37.genome.fa.gz -O /home/bioinformatikai/HW1/references/mm9.fa.gz

#Download the reference transcriptome in FASTA format:
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M1/gencode.vM1.pc_transcripts.fa.gz -O /home/bioinformatikai/HW1/references/mm9_RNA.fa.gz

#Download GTF/GFF3 file for your reference genome:
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz -O /home/bioinformatikai/HW1/references/mm9.gtf.gz

#Download raw FASTQ files:
prefetch -O /home/bioinformatikai/HW1/inputs SRR8985047 SRR8985048 SRR8985051 SRR8985052

fastq-dump --outdir /home/bioinformatikai/HW1/inputs/ --gzip --split-files /home/bioinformatikai/HW1/inputs/SRR8985047 /home/bioinformatikai/HW1/inputs/SRR8985048 /home/bioinformatikai/HW1/inputs/SRR8985051 /home/bioinformatikai/HW1/inputs/SRR8985052


