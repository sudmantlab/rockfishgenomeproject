#!/usr/bin/bash

set -e

# Basic bash for combining & converting the sequencing data for assembly later 

zcat Seb111_Cell*/*.subreads.fastq.gz >aleutianus_pacbio.fastq; gzip < aleutianus_pacbio.fastq >aleutianus_pacbio.fastq.gz

module load smrtlink/v600
dataset merge $PWD/aleutianus_pacbio.xml $PWD/*/*subreadset.xml

mkdir -p $PWD/split_fasta; 
python2.7 ~/.local/bin/FALCON-formatter $PWD/aleutianus_pacbio.fastq -o $PWD/split_fasta/

\ls -d1 $PWD/split_fasta/*.fasta >input_fasta.fofn; 
\ls -d1 $PWD/S*_Cell*/*subreads.bam >input_bam.fofn; 

