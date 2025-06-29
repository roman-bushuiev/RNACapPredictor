## Workflow

1. Define possible isoforms. Create an isoforms BLAST database (using distinct local regions).

```shell
makeblastdb -in references.fasta -dbtype nucl -out references_db
```

2. For each read in each `ba{X}.fasta` file, BLAST the most similar isoform. The output will be the file assigning the closest isoform for each read.

```shell
blastn -query bp01.fasta -db isoform_database/isoforms_reviewed_core_motifs -out blast_results/bp01_results.txt -outfmt "6 qseqid sseqid pident length evalue bitscore" -max_target_seqs 1 -num_threads 4
...
```

3. Run `BLAST_results_parsing.ipynb` to group reads by isoforms into separate files.

4. Convert `.lst` files (read IDs) from the previous step into `.fasta` files (read sequences).

```shell
seqtk subseq bp01.fasta separated_reads_lists/bp01_isoform_A_\(1\)_read_IDs.lst >       reads_separated/bp01_reads_from_isoform_A.fasta
...
```

5. Prepare isoforms for alignment with Bowtie2 by building Bowtie2 index for each isoform.

```shell
bowtie2-build ./individual_isoforms/isoform_A.fasta ./bowtie2/index/fm165_isoform_A
```

6. Align reads in `.fasta` files from step 4 using Bowtie2. The output is a `.bam` file for each `.fasta` file.

```
bowtie2 --local -f -p 8 -N 1 -L 6 -x ./bowtie2/index/fm165_isoform_A -U ./reads_separated/bp01_reads_from_isoform_A.fasta | samtools view -Shb | samtools sort -o bowtie2/FM165_isoform_A_bp01.bam
``

## Installation

conda install bioconda::blast
conda install bioconda::seqtk
conda install bioconda::bowtie2
conda install -c bioconda samtools