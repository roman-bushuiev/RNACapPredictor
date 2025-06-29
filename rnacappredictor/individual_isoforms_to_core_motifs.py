import os
from pathlib import Path
import re

def extract_core_motif(sequence):
    # Find pattern starting with constant sequence, followed by T's, followed by rest until N
    match = re.match(r'^CTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG[T]+(.+?)N', sequence)
    if match:
        return match.group(1)
    return None

def process_fasta_folder(input_folder):
    input_path = Path(input_folder)
    output_file = input_path.parent / 'isoforms_reviewed_core_motifs.fasta'
    
    with open(output_file, 'w') as outf:
        # Process each .fasta file in the folder
        for fasta_file in input_path.glob('*.fasta'):
            with open(fasta_file) as f:
                # Read header and sequence
                header = f.readline().strip()
                sequence = f.readline().strip()
                
                # Extract core motif
                core_motif = extract_core_motif(sequence)
                
                if core_motif:
                    # Write to output file
                    outf.write(f'{header}\n')
                    outf.write(f'{core_motif}\n')

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2:
        print("Usage: python individual_isoforms_to_core_motifs.py <input_folder>")
        sys.exit(1)
        
    input_folder = sys.argv[1]
    process_fasta_folder(input_folder)
