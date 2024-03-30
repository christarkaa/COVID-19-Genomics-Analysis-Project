#!/usr/bin/env bash

# Error handling: Exit script if any command fails
set -e

# Create a directory for the project
mkdir -p covid && cd covid || exit 1
echo "Created project directory."

# Download the dataset (assuming ERR5743893 is a valid accession)
fastq-dump --split-files ERR5743893
echo "Downloaded dataset."

# Create directories for quality control and mapping
mkdir -p QC_Reports
mkdir -p Mapping
echo "Created directories for quality control and mapping."

# Quality control
fastqc ERR5743893_1.fastq ERR5743893_2.fastq -o QC_Reports
echo "Performed quality control."

# Summarize QC results
multiqc QC_Reports
echo "Summarized QC results."

# Index the reference genome for mapping (if not already indexed)
if [ ! -f MN908947.fasta.bwt ]; then
    bwa index MN908947.fasta
    echo "Indexed reference genome for mapping."
fi

# Mapping
bwa mem MN908947.fasta ERR5743893_1.fastq ERR5743893_2.fastq > Mapping/ERR5743893.sam
samtools view -@ 20 -S -b Mapping/ERR5743893.sam > Mapping/ERR5743893.bam
samtools sort -@ 32 -o Mapping/ERR5743893.sorted.bam Mapping/ERR5743893.bam
samtools index Mapping/ERR5743893.sorted.bam
echo "Performed mapping."

# Index the reference genome for variant calling (if not already indexed)
if [ ! -f MN908947.fasta.fai ]; then
    samtools faidx MN908947.fasta
    echo "Indexed reference genome for variant calling."
fi

# Variant calling using freebayes
freebayes -f MN908947.fasta Mapping/ERR5743893.sorted.bam > ERR5743893.vcf
echo "Performed variant calling."

# Compress and index the VCF file
bgzip ERR5743893.vcf
tabix ERR5743893.vcf.gz
echo "Compressed and indexed VCF file."

echo "Analysis completed successfully!"




