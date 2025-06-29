#!/bin/bash

# SOURCE_DIR="data/FM179_allBPs/no_sample_id/20250204_1454_MN43023_FBA89620_52924708/fastq_pass"
# SOURCE_DIR="data/FM179_BP13-20/20250211_1617_MN43023_FAZ12209_ac671a55/fastq_pass"
# SOURCE_DIR="data/FM179_BP01-12/no_sample_id/20250205_1648_MN43023_FBA88286_f95b3399/fastq_pass"
# TARGET_DIR="data/FM179_combined/fastq_pass"

# # Loop over each barcode* directory in source
# for src_barcode_dir in "$SOURCE_DIR"/barcode*; do
#     # Get just the folder name (e.g., barcode01)
#     barcode_name=$(basename "$src_barcode_dir")
    
#     # Construct the corresponding target path
#     target_barcode_dir="$TARGET_DIR/$barcode_name"
    
#     # Create target barcode folder if it doesn't exist
#     mkdir -p "$target_barcode_dir"
    
#     # Move all files from source barcode folder to target
#     echo "Moving reads from $src_barcode_dir to $target_barcode_dir"
#     mv "$src_barcode_dir"/* "$target_barcode_dir"/
# done

#!/bin/bash

SOURCE_DIR="data/FM183/no_sample_id/20250306_1601_MN43023_FBC12602_17f7d59a/fastq_pass"
TARGET_DIR="data/FM183_combined_bps"

# Function to pad number with leading zero
pad_number() {
    printf "%02d" $1
}

# Process groups of 3 barcodes
for i in {1..5}; do
    # Calculate the three barcode numbers for this group
    b1=$(pad_number $i)
    b2=$(pad_number $((i+5)))
    b3=$(pad_number $((i+10)))
    
    # Create merged directory name
    merged_dir="$TARGET_DIR/barcode${b1}"
    mkdir -p "$merged_dir"
    
    # Move files from each barcode directory to merged directory
    for barcode in "barcode${b1}" "barcode${b2}" "barcode${b3}"; do
        if [ -d "$SOURCE_DIR/$barcode" ]; then
            echo "Moving files from $SOURCE_DIR/$barcode to $merged_dir"
            cp "$SOURCE_DIR/$barcode"/* "$merged_dir"/ 2>/dev/null || true
        fi
    done
done
