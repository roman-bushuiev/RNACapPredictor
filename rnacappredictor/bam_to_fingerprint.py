import pysam
from collections import Counter
import argparse
import os
import pandas as pd


def normalize_fingerprints_df(df):
    df['num_reads_ACGT'] = df['num_A'] + df['num_C'] + df['num_G'] + df['num_T']
    for l in ['A', 'C', 'G', 'T', 'INS', 'DEL']:
        df[f'{l}%_INSDEL'] = df[f'num_{l}'] / df['num_reads']
        df = df.fillna(0)
    for l in ['A', 'C', 'G', 'T']:
        df[f'{l}%'] = df[f'num_{l}'] / df['num_reads_ACGT']
        df = df.fillna(0)
    return df


def postprocess_fingerprints_df(df):
    print(df['file_name'])
    df.insert(0, 'barcode', df['file_name'].str.extract(r'(barcode\d+)')[0])
    df.insert(1, 'isoform', df['file_name'].str.extract(r'isoform_(.+)(?:\.bam)')[0])
    df = df.drop(columns=['file_name'])
    print(df['barcode'])
    df['barcode_num'] = df['barcode'].apply(lambda x: int(x.split('barcode')[1]))
    df = normalize_fingerprints_df(df)
    df = df.sort_values(by='barcode_num')
    return df


def analyze_position(bamfile, ref_name, pos):
    base_counts = Counter()
    strand_counts = {'+': Counter(), '-': Counter()}
    insertion_counts = Counter()
    total_reads = 0
    deletion_count = 0

    for read in bamfile.fetch(ref_name, pos, pos + 1):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        strand = '-' if read.is_reverse else '+'
        ref_pos = read.reference_start
        query_pos = 0

        for cigartype, length in read.cigartuples:
            if cigartype == 0:  # MATCH or MISMATCH (M)
                for i in range(length):
                    if ref_pos == pos:
                        base = read.query_sequence[query_pos].upper()
                        base_counts[base] += 1
                        strand_counts[strand][base] += 1
                        total_reads += 1
                    ref_pos += 1
                    query_pos += 1
            elif cigartype == 1:  # INSERTION (I)
                if ref_pos == pos:
                    ins_seq = read.query_sequence[query_pos:query_pos + length].upper()
                    insertion_counts[ins_seq] += 1
                    total_reads += 1
                query_pos += length
            elif cigartype == 2:  # DELETION (D)
                if ref_pos == pos:
                    deletion_count += 1
                    total_reads += 1
                ref_pos += length
            elif cigartype in [3, 4, 5]:  # skip: REF_SKIP, SOFT_CLIP, HARD_CLIP
                if cigartype in [4, 5]:
                    query_pos += length
                if cigartype == 3:
                    ref_pos += length

    return base_counts, strand_counts, insertion_counts, deletion_count, total_reads


def process_bam_file(bam_path, ref_path):
    """Process a single BAM file and return its statistics"""
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    ref = pysam.FastaFile(ref_path)

    ref_name = ref.references[0]
    ref_seq = ref.fetch(ref_name)

    n_positions = [i for i, base in enumerate(ref_seq) if base.upper() == 'N']
    if len(n_positions) != 1:
        print(f"Error: Expected 1 N position, but found {len(n_positions)} in {ref_name}")
        return None

    pos = n_positions[0]
    base_counts, _, insertions, deletions, total_reads = analyze_position(bamfile, ref_name, pos)

    bamfile.close()
    ref.close()

    return {
        'file_name': os.path.basename(bam_path),
        'num_reads': total_reads,
        'num_A': base_counts['A'],
        'num_C': base_counts['C'],
        'num_G': base_counts['G'],
        'num_T': base_counts['T'],
        'num_DEL': deletions,
        'num_INS': sum(insertions.values())
    }


def main():
    parser = argparse.ArgumentParser(description='Analyze BAM files for nucleotide distributions at N position')
    parser.add_argument('--data_folder', required=True, help='Folder containing BAM files')
    parser.add_argument('--out_pth', required=True, help='Path to save the output CSV file')
    parser.add_argument('--isoforms_dir', required=True, help='Path to the directory containing isoform fasta files')
    args = parser.parse_args()

    results = []
    
    # Process all BAM files in the data folder
    for filename in os.listdir(args.data_folder):
        if filename.endswith('.bam'):
            print(f"Processing {filename}")
            bam_path = os.path.join(args.data_folder, filename)
            # Construct reference path based on isoform name in BAM filename
            isoform = filename.split('_')[2].split('.')[0]  # Assuming format like "barcode01_isoform_A_reads.bam"
            ref_path = os.path.join(args.isoforms_dir, f'isoform_{isoform}.fasta')
            
            stats = process_bam_file(bam_path, ref_path)
            if stats:
                results.append(stats)
            print(stats)

    # Create DataFrame and save to CSV
    df = pd.DataFrame(results)
    df = postprocess_fingerprints_df(df)
    df.to_csv(args.out_pth, index=False)
    print(f"Results saved to {args.out_pth}")


if __name__ == "__main__":
    main()