#!/usr/bin/env bash

# Make a new directory for the project and change to the covid directory
mkdir covid && cd "$_"

# Download the dataset 
fastq-dump --split-files ERR5743893

# Make a new directory for quality control reports
mkdir -p QC_Reports 

# Quality control
fastqc ERR5743893_1.fastq ERR5743893_2.fastq -o QC_Reports

# Use MultiQC to summarise the QC results
multiqc . 

# The reads are good so there is no need for trimming

# Make a new directory "Mapping"
mkdir Mapping

# Download the reference genome
wget "https://www.futurelearn.com/links/f/no858mqqw7cxpdgv3ko0eqohoo2p7qc"

# Rename the file
mv MN908947.fasta.sa MN908947.fasta

# Index the reference genome
bwa index MN908947.fasta 

# Mapping
bwa mem MN908947.fasta ERR5743893_1.fastq ERR5743893_2.fastq > Mapping/ERR5743893.sam

# Convert the SAM file to a BAM file to save space
samtools view -@ 20 -S -b Mapping/ERR5743893.sam > Mapping/ERR5743893.bam

# Sort the BAM file based on the order the reads were mapped
samtools sort -@ 32 -o Mapping/ERR5743893.sorted.bam Mapping/ERR5743893.bam

# Index the BAM file
samtools index Mapping/ERR5743893.sorted.bam

# Index the reference genome using samtools
samtools faidx MN908947.fasta

# Variant calling using freebayes
freebayes -f MN908947.fasta Mapping/ERR5743893.sorted.bam  > ERR5743893.vcf

# Compres the VCF file
bgzip ERR5743893.vcf

# Index the VCF file
tabix ERR5743893.vcf.gz




