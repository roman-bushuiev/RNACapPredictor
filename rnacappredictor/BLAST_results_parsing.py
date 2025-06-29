def main():
    import argparse
    import statistics
    from statistics import StatisticsError
    import glob
    import os

    # Set up command line arguments
    parser = argparse.ArgumentParser(description='Parse BLAST results and create files with readnames by isoform')
    parser.add_argument('--blast_dir', type=str, required=True, help='Path to directory containing BLAST results')
    parser.add_argument('--isoforms_dir', type=str, required=True, help='Path to directory containing isoforms')
    args = parser.parse_args()

    # Get all BLAST result files
    blast_files = sorted(glob.glob(os.path.join(args.blast_dir, '*_results.txt')))

    #script parses best hit BLAST results and creates files with readnames belonging to a given isoform
    for blast_file in blast_files:
        # Extract barcode/sample ID from filename
        current_bp = os.path.basename(blast_file).split('_')[0]

        # each element has [isoform-name, [read-names], [bit-scores], [e-value]]
        ## bitscore is initiated as 100, this is not correct, but in the grand scheem of things means very small error
        # Initialize empty dict
        isoform_dict = {}
        
        # Get all fasta files from isoforms directory and create dict entries
        for fasta_file in glob.glob(os.path.join(args.isoforms_dir, '*.fasta')):
            isoform_name = os.path.splitext(os.path.basename(fasta_file))[0]
            isoform_dict[isoform_name] = [[], [], []]

        with open(blast_file, 'r') as inFile:
            for line in inFile:
                readID = line.split()[0]
                isoformID = line.split()[1]
                readBIT = float(line.split()[5])
                readEVAL = float(line.split()[4])
                
                #add the readID
                isoform_dict[isoformID][0].append(readID)
                
                #add the BIT score
                isoform_dict[isoformID][1].append(readBIT)
                
                #add the E-value
                isoform_dict[isoformID][2].append(readEVAL)

        print(f"\nResults for {current_bp}:")
        for key in isoform_dict.keys():
            
            #BIT average
            try:
                bit_average = round(statistics.mean(isoform_dict[key][1]),2)
            except StatisticsError:
                bit_average = 'NA'
            
            #EVAL average
            try:
                eval_average = statistics.mean(isoform_dict[key][2])
            except StatisticsError:
                eval_average = 'NA'
            
            print('{:<20}\t{:<7}\t\t{:<7}\t\t{:<10}\n'.format(key,
                                                        len(isoform_dict[key][0]),
                                                        bit_average,
                                                        eval_average))

        # Save names of reads for each isoform found in isoform_dict
        for isoform_name in isoform_dict.keys():
            output_path = os.path.join(args.blast_dir, f'{current_bp}_{isoform_name}_read_IDs.lst')
            with open(output_path, 'w') as outFile:
                for readID in isoform_dict[isoform_name][0]:
                    outFile.write(f'{readID}\n')


if __name__ == '__main__':
    main()
