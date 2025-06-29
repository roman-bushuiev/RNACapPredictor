#! /bin/bash

if [ -z "$data_folder" ]; then
    echo "Error: data_folder environment variable is not set"
    exit 1
fi

if [ -z "$bowtie2_isoforms_folder" ]; then
    echo "Error: bowtie2_isoforms_folder environment variable is not set"
    exit 1
fi

echo "data_folder: $data_folder"
echo "bowtie2_isoforms_folder: $bowtie2_isoforms_folder"

blast_results_folder="$data_folder/blast_results"
bowtie2_results_folder="$data_folder/bowtie2_results"
mkdir -p "$bowtie2_results_folder"

for reads_fasta in "$blast_results_folder"/*_reads.fasta; do
    if [ -f "$reads_fasta" ]; then
        echo "Processing $reads_fasta"

        # Extract filename components
        filename=$(basename "$reads_fasta" _reads.fasta)
        barcode=$(echo "$filename" | cut -d'_' -f1)
        isoform=$(echo "$filename" | grep -o "isoform_[A-Z0-9-]*")
        
        echo "--------------------------------"
        echo "Fasta file: $reads_fasta"
        echo "Filename: $filename"
        echo "Barcode: $barcode"
        echo "Isoform: $isoform"

        # Align reads to the isoform
        # -f: output in FASTA format
        # -p 8: use 8 threads
        # -N 1: maximum number of alignments to report
        # -L 6: minimum alignment length
        bowtie2 \
            --local \
            -f \
            -p 8 \
            -N 1 \
            -L 6 \
            -x "$bowtie2_isoforms_folder/$isoform" \
            -U "$reads_fasta" \
            | samtools view -Shb \
            | samtools sort -o "$bowtie2_results_folder/${filename}.bam"

        # Index the BAM file (generate .bam.bai file from .bam file)
        samtools index "$bowtie2_results_folder/${filename}.bam"
    fi
done
