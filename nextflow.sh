#!/usr/bin/env bash

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


