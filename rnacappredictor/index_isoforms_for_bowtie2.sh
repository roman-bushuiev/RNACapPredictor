#! /bin/bash

# export isoforms_dir="data/isoforms/FM200"

if [ -z "$isoforms_dir" ]; then
    echo "Error: isoforms_dir environment variable is not set"
    exit 1
fi

mkdir -p $isoforms_dir/bowtie2_index

for fasta_file in $isoforms_dir/individual_isoforms/*.fasta; do
    # Extract the isoform name without path and extension
    isoform_name=$(basename "$fasta_file" .fasta)

    # Build Bowtie2 index
    bowtie2-build "$fasta_file" "$isoforms_dir/bowtie2_index/$isoform_name"
done
