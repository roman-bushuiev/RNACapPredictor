#!/bin/bash

if [ -z "$data_folder" ]; then
    echo "Error: data_folder environment variable is not set"
    exit 1
fi

if [ -z "$isoform_database" ]; then
    echo "Error: isoform_database environment variable is not set"
    exit 1
fi

echo "data_folder: $data_folder"
echo "isoform_database: $isoform_database"

blast_results_folder="$data_folder/blast_results"
mkdir -p "$blast_results_folder"

# 1. Process each FASTA file in the data folder and create BLAST results .txt files
for fasta_file in "$data_folder"/*.fasta; do
    if [ -f "$fasta_file" ]; then
        # Get the filename without path and extension
        filename=$(basename "$fasta_file" .fasta)
        echo "BLASTing $filename..."
        
        # Create output file path
        output_file="$blast_results_folder/${filename}_results.txt"
        
        # Run BLAST
        blastn -query "$fasta_file" \
              -db "$isoform_database" \
              -out "$output_file" \
              -outfmt "6 qseqid sseqid pident length evalue bitscore" \
              -max_target_seqs 1 \
              -num_threads 4

        # Check if BLAST completed successfully
        if [ -f "$output_file" ]; then
            result_count=$(wc -l < "$output_file")
            echo "BLAST completed for $filename - $result_count alignments found"
        else
            echo "Error: BLAST failed for $filename"
        fi
    fi
done

# 2. Parse BLAST .txt results and create .lst files with read names for each isoform
python3 rnacappredictor/BLAST_results_parsing.py --blast_dir "$blast_results_folder" --isoforms_dir "$isoforms_dir"

# 3. Iterate over all .lst files and create corresponding .fasta files
for lst_file in "$blast_results_folder"/*_read_IDs.lst; do
    if [ -f "$lst_file" ]; then
        echo "Processing $lst_file"

        # Extract barcode from filename (e.g. barcode01 from barcode01_isoform_A_(1)_read_IDs.lst)
        barcode=$(echo $(basename "$lst_file") | cut -d'_' -f1)
        
        # Create output filename
        output_fasta="${lst_file%_read_IDs.lst}_reads.fasta"
        
        # Run seqtk subseq using the appropriate barcode fasta file
        seqtk subseq "$data_folder/${barcode}.fasta" "$lst_file" > "$output_fasta"
        
        echo "Created $output_fasta"
    fi
done
