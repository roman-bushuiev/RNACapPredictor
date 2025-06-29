#!/bin/bash
# Concatenate all fastq.gz files in each barcode subfolder of the fastq_pass_folder into a single fastq file

if [ -z "$data_folder" ]; then
    echo "Error: data_folder environment variable is not set"
    exit 1
fi

echo "data_folder: $data_folder"

cd $data_folder

# Only process directories that start with 'barcode'
for D in barcode*/; do
    if [ -d "$D" ]; then
        echo "Processing $D"
        cd "$D"
        # Count the number of files to concatenate
        file_count=$(ls *.fastq.gz | wc -l)
        echo "Found $file_count .fastq.gz files in $D"
        
        # Concatenate the files
        zcat *.fastq.gz > ../"${D%/}.fastq"
        
        # Verify the output and count reads
        if [ -f "../${D%/}.fastq" ]; then
            total_lines=$(wc -l < "../${D%/}.fastq")
            read_count=$((total_lines / 4))
            echo "Successfully created ../${D%/}.fastq with $read_count reads ($total_lines lines)"
            
            # Convert FASTQ to FASTA using seqtk
            seqtk seq -a "../${D%/}.fastq" > "../${D%/}.fasta"
            
            # Verify FASTA file was created
            if [ -f "../${D%/}.fasta" ]; then
                fasta_lines=$(wc -l < "../${D%/}.fasta")
                fasta_reads=$((fasta_lines / 2))
                echo "FASTA file created with $fasta_reads reads ($fasta_lines lines)"
            else
                echo "Error: Failed to create FASTA file"
            fi
        else
            echo "Error: Failed to create ../${D%/}.fastq"
        fi
        cd ..
    fi
done 