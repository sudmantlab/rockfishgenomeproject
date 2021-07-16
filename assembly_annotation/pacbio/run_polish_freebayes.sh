#!/usr/bin/bash

set -e

##### Modified for local runs from VGP_assembly github - https://github.com/VGP/vgp-assembly

if [ $# -ne 3 ]; then
	printf "\n\n\tPlease provide 3 parameters as follows\n";
	printf "\t\trun_raconpolishillumina.sh FASTQ1 FASTQ2 REF_FASTA\n\n";
	exit 1;
fi

fastq1=$1
fastq2=$2
fasta=$3

reference=$(echo $fasta | sed -e 's/.*\///')

module load minimap2 samtools java

bam=$(echo $reference | sed -e 's/.*\///' -e 's/\.fasta/.1.bam/' -e 's/\.fas/.1.bam/' -e 's/\.fa/.1.bam/')
bam=$(echo $reference | sed -e 's/.*\///' -e 's/\.fasta/.2.bam/' -e 's/\.fas/.2.bam/' -e 's/\.fa/.2.bam/')
consensus=$(echo $reference | sed -e 's/.*\///' -e 's/\.fa$/.frby/' -e 's/\.fas$/.frby/' -e 's/\.fasta$/.frby/')

printf "\tMapping short reads\n";
minimap2 -t 24 -ax sr ${fasta} ${fastq1} ${fastq2} | samtools view -bh >${reference}.bam;
printf "\tSorting mapped reads\n";
samtools sort -@ 24 -T ./ -o ${reference}.sort.bam ${reference}.bam; rm ${reference}.bam &

printf "\tRunning GATK for duplicates\n";
Software/gatk-4.1.2/gatk --java-options "-Xmx100G" MarkDuplicates --TMP_DIR $PWD --INPUT ${reference}.sort.bam --OUTPUT ${bam} --METRICS_FILE ${bam}.metrics --MAX_FILE_HANDLES 15000;
rm ${reference}.sort.bam &

module load samtools freebayes vcftools gcc bcftools
samtools index ${bam} && samtools faidx ${fasta};
printf "\tFreebayesing now\n";
freebayes --max-coverage 200 -f ${fasta} ${bam} | bcftools view -Ou -e'type="ref"' | bcftools norm -Ob -f ${fasta} -o ${bam/%.bam/.bcf} --threads 32;
bcftools index ${bam/%.bam/.bcf};
bcftools consensus -i 'QUAL>1 && (GT="AA" || GT="Aa")' -Hla -f $fasta ${bam/%.bam/.bcf} > ${consensus}.1.fasta 
bcftools view -i 'QUAL>1 && (GT="AA" || GT="Aa")' -Oz --threads=32 ${bam/%.bam/.bcf} > ${bam/%.bam/.vcf.gz}

printf "\t\t\nROUND 2\n\n\tMapping short reads\n";
minimap2 -t 24 -ax sr ${consensus}.1.fasta ${fastq1} ${fastq2} | samtools view -bh >${reference}.bam;
printf "\tSorting mapped reads\n";
samtools sort -@ 24 -T ./ -o ${reference}.sort2.bam ${reference}.bam; rm ${reference}.bam &

printf "\tRunning GATK for duplicates\n";
Software/gatk-4.1.2/gatk --java-options "-Xmx100G" MarkDuplicates --TMP_DIR $PWD --INPUT ${reference}.sort2.bam --OUTPUT ${bam2} --METRICS_FILE ${bam2}.metrics --MAX_FILE_HANDLES 15000;
rm ${reference}.sort2.bam &

module load samtools freebayes vcftools gcc bcftools
samtools index ${bam2} && samtools faidx ${consensus}.1.fasta;
printf "\tFreebayesing now\n";
freebayes --max-coverage 200 -f ${consensus}.1.fasta ${bam2} | bcftools view -Ou -e'type="ref"' | bcftools norm -Ob -f ${consensus}.1.fasta -o ${bam2/%.bam/.bcf} --threads 32;
bcftools index ${bam2/%.bam/.bcf};
bcftools consensus -i 'QUAL>1 && (GT="AA" || GT="Aa")' -Hla -f ${consensus}.1.fasta ${bam2/%.bam/.bcf} > ${consensus}.2.fasta
bcftools view -i 'QUAL>1 && (GT="AA" || GT="Aa")' -Oz --threads=32 ${bam2/%.bam/.bcf} > ${bam2/%.bam/.vcf.gz}


