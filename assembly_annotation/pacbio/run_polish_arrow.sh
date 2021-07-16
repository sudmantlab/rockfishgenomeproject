#!/usr/bin/bash

set -e

# Using ARROW from smrtlink for polishing using pacbio reads

if [ $# -ne 3 ]; then
	printf "\n\n\tPlease provide 3 parameters as follows\n";
	printf "\t\trun_arrow.sh FASTQ REF_FASTA BAM_FOFN\n\n";
	exit 1;
fi

fastq=$1
fastaseq=$2
bam_fofn=$3

reference=$(echo $fastaseq | sed -e 's/.*\///' -e 's/\.fasta$//' -e 's/\.fas$//' -e 's/\.fa$//' -e 's/\.FASTA$//' -e 's/\.FAS$//' -e 's/\.FA$//')

conda activate pb-assembly
module load minimap2 samtools #smrtlink

count=0
ln -sf $fastaseq ${reference}.${count}.fasta
samtools faidx ${reference}.${count}.fasta

uniqID=$(echo $BASHPID)

for iter in 1; do
    if [ -e "${reference}.t1.fasta" ]; then
        printf "\t${reference}.t1.fasta EXISTS, hence skipping arrow polishing\n";
        exit 0;
    fi    
	if [ ! -e "pbmap_${reference}.${count}.bam" ]; then
		pbmm2 align --preset SUBREAD ${reference}.${count}.fasta $fastq /clusterfs/genomicdata/rockfish/tmp/pbmap_${reference}.${count}.bam --sort -j 26 -J 6;
        mv /clusterfs/genomicdata/rockfish/tmp/pbmap_${reference}.${count}.bam .;
        printf "\t pbmap_${reference}.${count}.bam mapped\n";
    fi
    if [ ! -e "pbmap_${reference}.pb.bam" ]; then
        printf "\t pbmap_${reference}.${count}.bam conversion to arrow compatible\n";
        pbbamify --input=pbmap_${reference}.${count}.bam --output=pbmap_${reference}.pb.bam ${reference}.${count}.fasta $bam_fofn
        printf "\t pbmap_${reference}.pb.bam modified\n";
	fi	
	if [ ! -e "pbmap_${reference}.pbsort.bam" ]; then
        printf "\t pbmap_${reference}.pbsort.bam sorting\n";
        samtools sort -T /clusterfs/genomicdata/rockfish/tmp/${uniqID} -@ 20 -m 1G -o pbmap_${reference}.pbsort.bam pbmap_${reference}.pb.bam
		pbindex pbmap_${reference}.pbsort.bam;
        printf "\t pbmap_${reference}.pbsort.bam INDEXED\n";
	fi	
	printf "\tVariant calling pbmap_${reference}.pbsort.bam\n";
    variantCaller pbmap_${reference}.pbsort.bam --algorithm=arrow -j 20 -r ${reference}.${count}.fasta -o ${reference}.t1.fasta;
	samtools faidx ${reference}.t1.fasta;
	count=$((count+1));
done

printf "\tPolished to produce ${reference}.t1.fasta\n";

