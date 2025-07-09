#! /bin/bash

# export isoforms_dir="data/isoforms/FM200"

if [ -z "$isoforms_dir" ]; then
    echo "Error: isoforms_dir environment variable is not set"
    exit 1
fi

makeblastdb -in $isoforms_dir/isoforms_reviewed_core_motifs.fasta -dbtype nucl -out $isoforms_dir/isoform_database/isoform_database