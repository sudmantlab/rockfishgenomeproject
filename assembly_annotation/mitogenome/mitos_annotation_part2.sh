#!/usr/bin/bash

set -e

\ls -d /global/scratch2/Rockfish/Genomes/sequencing/illumina/output/*/*/mitogenome/MITOS/sequence.fas-8* | grep -v 'work' |
    while read fasta2; do

        working=$(echo $fasta2 | sed -e 's/\/mitogenome\/.*/\/mitogenome\//');
        speciesname=$(echo $fasta2 | sed -e 's/\/mitogenome\/.*//' -e 's/.*\/output\///' -e 's/\/.*//');
        samplename=$(echo $fasta2 | sed -e 's/\/mitogenome\/.*//' -e 's/.*\/output\///' -e 's/.*\///');
        name=$(echo $working | sed -e 's/\/mitogenome\/.*//' -e 's/.*\/output\///' -e 's/\//./');
        fasta1=$(\ls -d ${working}/*.result/*.fasta | grep -v -e 'work' -e '\.fai');

        printf "Running MITOS2\n";
        printf "Working at ${working} \n \tOn ${fasta1}\n";
        module load bedtools samtools;
        samtools faidx ${fasta1};

        fgrep -w -e 'trnF(gaa)' ${working}/MITOS/result.bed >${working}/MITOS/phe.tmp;
        nonproper=$(awk -F'\t' 'NR==1 {if($NF=="-") {print "YES"} else print "NO"}' ${working}/MITOS/phe.tmp);     ### STRAND REVERSAL
        seqname=$(awk '{print $1}' ${working}/MITOS/phe.tmp);
        phe_start=$(awk '{print $2-1}' ${working}/MITOS/phe.tmp);
        phe_end=$(awk '{print $3}' ${working}/MITOS/phe.tmp);
        size_pre=$(awk '{print $3}' ${working}/MITOS/phe.tmp);
        size_end=$(awk 'NR==1 {print $2}' ${fasta1}.fai);
        rm ${working}/MITOS/phe.tmp;

        printf "${seqname}\t${phe_start}\t${phe_end}\n${seqname}\t${size_pre}\t${size_end}\n${seqname}\t0\t${phe_start}\n" \
            >${working}/ext.bed ;                           ### Creating coordinates for extracting based on trnPhe

        printf ">${speciesname}___${samplename}\n" \
            >${working}/${name}.fasta;
        rm -fr ${working}/MITOS2/*;

        if [ "${nonproper}" == "YES" ]; then                ### Negative strand for trn-Phe
            printf "\tRunning in Antisense direction on\n${fasta1}\n";
            tac ${working}/ext.bed | 
                bedtools getfasta -fi ${fasta1} -bed - |
                grep -v '^>' | 
                tr "\n" " " |
                sed -e 's/ //g' | 
                tr "[ATGCatgc]" "[TACGtacg]" | rev |
                sed -e 's/$/\n/' \
                    >>${working}/${name}.fasta;
            mkdir -p ${working}/MITOS2;
            printf "\tRunning MITOS local version\n";
            runmitos.py -i ${working}/${name}.fasta -o ${working}/MITOS2 \
                --refdir /global/scratch2/miniconda3/envs/mitos/lib/python2.7/site-packages/mitos/data/refseq81m \
                -c 2 --evalue 1e-3 --best --ncev 1e-2;

        else
            printf "\tRunning in Sense direction on\n${fasta1}\n";
            bedtools getfasta -fi ${fasta1} -bed ${working}/ext.bed |
                grep -v '^>' |
                tr "\n" " " |
                sed -e 's/ //g' -e 's/$/\n/' \
                    >>${working}/${name}.fasta;
            mkdir -p ${working}/MITOS2;
            printf "\tRunning MITOS local version\n";
            runmitos.py -i ${working}/${name}.fasta -o ${working}/MITOS2 \
                --refdir /global/scratch2/miniconda3/envs/mitos/lib/python2.7/site-packages/mitos/data/refseq81m \
                -c 2 --evalue 1e-3 --best --ncev 1e-2;        
        fi   
        samtools faidx ${working}/${name}.fasta;
   done    


