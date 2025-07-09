#!/bin/bash

export isoforms_dir="data/isoforms/U1-11"

python3 rnacappredictor/individual_isoforms_to_core_motifs.py $isoforms_dir/individual_isoforms
./rnacappredictor/index_isoforms_for_bowtie2.sh
./rnacappredictor/index_isoforms_for_BLAST.sh