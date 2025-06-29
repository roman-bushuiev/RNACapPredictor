#! /bin/bash

if [ -z "$data_folder" ]; then
    echo "Error: data_folder environment variable is not set"
    exit 1
fi

if [ -z "$isoforms_dir" ]; then
    echo "Error: isoforms_dir environment variable is not set" 
    exit 1
fi

echo "data_folder: $data_folder"
echo "isoforms_dir: $isoforms_dir"

python3 rnacappredictor/bam_to_fingerprint.py \
    --data_folder "$data_folder/bowtie2_results" \
    --out_pth "$data_folder/fingerprints.csv" \
    --isoforms_dir "$isoforms_dir"

