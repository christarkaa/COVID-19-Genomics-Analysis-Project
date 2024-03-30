#!/usr/bin/env bash

# Define an array of sample accession numbers
accessions=(
  "ERR5743893"
  "ERR5556343"
  "SRR13500958"
  "ERR5181310"
  "ERR5405022"
)

# Loop through each accession number
for accession in "${accessions[@]}"; do
  # Make a new directory for the project and change to the directory
  mkdir -p "$accession" && cd "$accession"

  # Download the dataset 
  fastq-dump --split-files "$accession"

  # Make a new directory for quality control reports
  mkdir -p QC_Reports 

  # Quality control
  fastqc "${accession}_1.fastq" "${accession}_2.fastq" -o QC_Reports

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
  bwa mem MN908947.fasta "${accession}_1.fastq" "${accession}_2.fastq" > "Mapping/${accession}.sam"

  # Convert the SAM file to a BAM file to save space
  samtools view -@ 20 -S -b "Mapping/${accession}.sam" > "Mapping/${accession}.bam"

  # Sort the BAM file based on the order the reads were mapped
  samtools sort -@ 32 -o "Mapping/${accession}.sorted.bam" "Mapping/${accession}.bam"

  # Index the BAM file
  samtools index "Mapping/${accession}.sorted.bam"

  # Index the reference genome using samtools
  samtools faidx MN908947.fasta

  # Variant calling using freebayes
  freebayes -f MN908947.fasta "Mapping/${accession}.sorted.bam"  > "${accession}.vcf"

  # Compress the VCF file
  bgzip "${accession}.vcf"

  # Index the VCF file
  tabix "${accession}.vcf.gz"

  # Convert the VCF file to a .csv file for further analysis and visualization
  echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER" > "${accession}.csv"
  bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\n' "${accession}.vcf.gz" >> "${accession}.csv"

  # Move back to the parent directory
  cd ..
done

