#!/usr/bin/bash

#####	02/04/2019	#####
#####	Conversion of subreads bam file to FASTQ	#####

set -e 

if [ "$#" -ne 1 ]; then
    printf "Illegal number of parameters:\t	$#	\n";
    printf "\tPlease RUN as : sh extract_subreads_bam2fastq.sh FilePattern\n";
    printf "\t\tExample:	sh extract_subreads_bam2fastq.sh /PATH/TO/SEB10_Cell\n";
    exit;
fi

SUBREADS_PATHING=$1"*/*.subreads.bam"

checking=$(ls $SUBREADS_PATHING | wc -l)
if [ "$checking" -eq 0 ]; then
    printf ". files with pattern $checking dont exist\n\tPLEASE CHECK for $1*/*.subreads.bam\n";
    exit
fi

for file1 in $SUBREADS_PATHING; do
  file2=$(echo $file1 | sed -e 's/\.bam$//'); 
  printf "Starting $file1\n";
  file_fq=$(echo $file1 | sed -e 's/\.bam$/.fq.gz/');

  depth=$(ls -lh $file1 | cat | awk '{print $NF}' | sed -e 's/\//\t/g' | awk -F'\t' '{print NF-1}' | sort -un | wc -l)

  if [ "$depth" -eq "0" ]; then
    printf "Something wrong with folder depth; Check if there are other subfolders with files in the following output\n";
    ls -lh $file1;
    exit;
  fi

  printf "  Running:\t Bam2Fastq for $file1\n\n";
  conversion=" /global/scratch2/Software/smrtlink/v600/smrtlink/smrtcmds/bin/bam2fastq -o $file2 $file1"
  printf " $conversion \n";
  eval $conversion


  printf "\tProcessed   $file1 \n"; 
done	


