#!/usr/bin/env bash

# Define an array of sample accession numbers
accessions=(
  "ERR5743893"
  "ERR5556343"
  "SRR13500958"
  "ERR5181310"
  "ERR5405022"
)

# Error handling: Exit script if any command fails
set -eux

# Define a function to process each sample
process_sample() {
    local accession="$1"

    # Create a directory for the sample
    mkdir -p "$accession"
    cd "$accession" || exit 1

    # Download the dataset 
    fastq-dump --split-files "$accession"

    # Create a directory for quality control reports
    mkdir -p QC_Reports 

    # Perform quality control
    fastqc "${accession}_1.fastq" "${accession}_2.fastq" -o QC_Reports

    # Use MultiQC to summarize the QC results
    multiqc . 

    # Create a directory for mapping
    mkdir -p Mapping

    # Download the reference genome if not already downloaded
    if [ ! -f MN908947.fasta ]; then
        wget "https://www.futurelearn.com/links/f/no858mqqw7cxpdgv3ko0eqohoo2p7qc"
        mv MN908947.fasta.sa MN908947.fasta
        bwa index MN908947.fasta 
        samtools faidx MN908947.fasta
    fi

    # Map the reads
    bwa mem MN908947.fasta "${accession}_1.fastq" "${accession}_2.fastq" > "Mapping/${accession}.sam"
    samtools view -@ 20 -S -b "Mapping/${accession}.sam" > "Mapping/${accession}.bam"
    samtools sort -@ 32 -o "Mapping/${accession}.sorted.bam" "Mapping/${accession}.bam"
    samtools index "Mapping/${accession}.sorted.bam"

    # Variant calling
    freebayes -f MN908947.fasta "Mapping/${accession}.sorted.bam" > "${accession}.vcf"
    bgzip "${accession}.vcf"
    tabix "${accession}.vcf.gz"

    # Convert VCF to CSV
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\n' "${accession}.vcf.gz" > "${accession}.csv"

    # Move back to the parent directory
    cd ..
}

# Process each sample
for accession in "${accessions[@]}"; do
    process_sample "$accession"
done
