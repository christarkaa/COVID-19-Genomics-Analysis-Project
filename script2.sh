#!/usr/bin/env bash

# Error handling: Exit script if any command fails
set -e

# Define an array of sample accession numbers
samples=(
  "ERR5743893"
  "ERR5556343"
  "SRR13500958"
  "ERR5181310"
  "ERR5405022"
)

# Index the reference genome for mapping (if not already indexed)
if [ ! -f MN908947.fasta.bwt ]; then
    bwa index MN908947.fasta
    echo "Indexed reference genome for mapping."
fi

# Loop through each sample
for sample in "${samples[@]}"; do
  echo "Processing sample: $sample"

  # Download the dataset (assuming ERR5743893 is a valid accession)
  fastq-dump --split-files "$sample"
  echo "Downloaded dataset for $sample."

  # Create directories for quality control and mapping
  mkdir -p "$sample"/QC_Reports
  mkdir -p "$sample"/Mapping
  echo "Created directories for quality control and mapping for $sample."

  # Quality control
  fastqc "${sample}"_1.fastq "${sample}"_2.fastq -o "$sample"/QC_Reports
  echo "Performed quality control for $sample."

  # Summarize QC results
  multiqc "$sample"/QC_Reports
  echo "Summarized QC results for $sample."

  # Mapping
  bwa mem MN908947.fasta "${sample}"_1.fastq "${sample}"_2.fastq > "$sample"/Mapping/"$sample".sam
  samtools view -@ 20 -S -b "$sample"/Mapping/"$sample".sam > "$sample"/Mapping/"$sample".bam
  samtools sort -@ 32 -o "$sample"/Mapping/"$sample".sorted.bam "$sample"/Mapping/"$sample".bam
  samtools index "$sample"/Mapping/"$sample".sorted.bam
  echo "Performed mapping for $sample."

  # Index the reference genome for variant calling (if not already indexed)
  if [ ! -f MN908947.fasta.fai ]; then
      samtools faidx MN908947.fasta
      echo "Indexed reference genome for variant calling."
  fi

  # Variant calling using freebayes
  freebayes -f MN908947.fasta "$sample"/Mapping/"$sample".sorted.bam > "$sample"/"$sample".vcf
  echo "Performed variant calling for $sample."

  # Compress and index the VCF file
  bgzip "$sample"/"$sample".vcf
  tabix "$sample"/"$sample".vcf.gz
  echo "Compressed and indexed VCF file for $sample."

  # Convert VCF to CSV
  echo -e "Sample\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER" > "$sample"/"$sample".csv
  bcftools query -l "$sample"/"$sample".vcf.gz | sed 's/$/'$'\t''/' > "$sample"/sample_column.txt
  paste "$sample"/sample_column.txt <(bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\n' "$sample"/"$sample".vcf.gz) >> "$sample"/"$sample".csv
  echo "Converted VCF to CSV for $sample."

  # Move the CSV file to the main directory
  mv "$sample"/"$sample".csv .

done

# Merge CSV files for all samples
echo "Merging CSV files..."
cat *.csv > merged.csv
echo "Merged CSV files for all samples."

echo "Analysis completed successfully!"
