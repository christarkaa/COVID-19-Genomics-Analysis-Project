#!/usr/bin/env bash

set -e

# Create a sample text file and copy and paste the Ascencion numbers
# Define the ENA accessions
accessions=(
  "ERR5556343"
  "SRR13500958"
  "ERR5743893"
  "ERR5181310"
  "ERR5405022"
)

# Create samples.txt and copy accessions into it
echo "${accessions[@]}" > samples.txt
echo "File samples.txt created with the ENA accessions."

# Use a for loop to run fastq-dump on each of the accessions in samples.txt
for i in $(cat samples.txt);do fastq-dump --split-files $i;done

# Compress the .fasta files
gzip *.fastq

# Create a directory and move the .fastq.gz files into it
mkdir data
mv *.fastq.gz data

#  Create the samplesheet file
# Download the python script for downloading creating the samplesheet file
wget -L https://raw.githubusercontent.com/nf-core/viralrecon/master/bin/fastq_dir_to_samplesheet.py

# Now, create the samplesheet file
python3 fastq_dir_to_samplesheet.py data samplesheet.csv -r1 _1.fastq.gz -r2 _2.fastq.gz

# view samplesheet file 
cat samplesheet.csv

# Run viralrecon
nextflow run nf-core/viralrecon -profile docker \
--max_memory '12.GB' --max_cpus 4 \
--input samplesheet.csv \
--outdir results/viralrecon \
--protocol amplicon \
--genome 'MN908947.3' \
--primer_set artic \
--primer_set_version 3 \
--skip_kraken2 \
--skip_assembly \
--skip_pangolin \
--skip_nextclade \
--platform illumina


