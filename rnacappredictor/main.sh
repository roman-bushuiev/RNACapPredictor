#! /bin/bash

# Example usage: ./rnacappredictor/main.sh 2>&1 | tee logs/FM165.log

# export data_folder="data/FM180_BP01-12/no_sample_id/20250224_1757_MN43023_FAZ12635_68f383d0/fastq_pass"
# export data_folder="data/FM181_BP13-20/no_sample_id/20250226_2018_MN43023_FAZ13458_a45f3d18/fastq_pass"
# export data_folder="data/FM179_combined/fastq_pass"
# export data_folder="data/FM185/no_sample_id/20250324_1822_MD-101425_FBC23952_d258aa5b/fastq_pass"
# export data_folder="data/FM165/no_sample_id/20241120_1928_MN43023_FAY46018_68f388f5/fastq_pass"
# export data_folder="data/FM183/no_sample_id/20250306_1601_MN43023_FBC12602_17f7d59a/fastq_pass"
# export data_folder="data/FM183_combined_bps"
# export data_folder="data/FM185/no_sample_id/20250324_1822_MD-101425_FBC23952_d258aa5b/fastq_pass"
# export data_folder="data/FM200/no_sample_id/20250609_1143_MD-101425_FBC21923_a6c9a0a2/fastq_pass"

export data_folder="data/FM200/no_sample_id/20250609_1143_MD-101425_FBC21923_a6c9a0a2/fastq_pass/U6"
export isoforms_name="U6"

export isoforms_dir="data/isoforms/$isoforms_name/individual_isoforms"
export isoform_database="data/isoforms/$isoforms_name/isoform_database/isoform_database"
export bowtie2_isoforms_folder="data/isoforms/$isoforms_name/bowtie2_index"

echo "Running 1_concat_raw_fastq_folders_into_fastas.sh"
./rnacappredictor/1_concat_raw_fastq_folders_into_fastas.sh || { echo "1_concat_raw_fastq_folders_into_fastas.sh failed"; exit 1; }

echo "Running 2_split_fastas_by_isoforms.sh"
./rnacappredictor/2_split_fastas_by_isoforms.sh || { echo "2_split_fastas_by_isoforms.sh failed"; exit 1; }

echo "Running 3_align_within_isoforms.sh"
./rnacappredictor/3_align_within_isoforms.sh || { echo "3_align_within_isoforms.sh failed"; exit 1; }

echo "Running 4_compute_fingerprints.sh"
./rnacappredictor/4_compute_fingerprints.sh || { echo "4_compute_fingerprints.sh failed"; exit 1; }

echo "Done"
