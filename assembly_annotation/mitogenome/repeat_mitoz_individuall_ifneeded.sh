#!/usr/bin/bash

set -e

maindir="Rockfish/Genomes/sequencing/illumina"

# We used error corrected Illumina PE data with Musket if needed (http://musket.sourceforge.net/homepage.htm)

cd ${maindir};

for fastq1 in ${maindir}/output/*_*/*/corrected/corr.0; do 

    fastq2=$(echo $fastq1 | sed -e 's/\/corr\.0/\/corr.1/');
    workdir=$(echo $fastq1 | sed -e 's/\/corrected\/corr.0$/\/mitogenome/');
    name=$(echo $fastq1 | sed -e 's/\/corrected\/corr.0$//' -e 's/.*\/output\///' -e 's/.*\///');
    resultdir=$(echo ${workdir} | sed -e 's/$/\/${name}.result/');
    mkdir -p ${workdir} &&  cd ${workdir};
    printf "Starting run for ${name} \n";
    if [ -e "${workdir}/${name}.result/${name}.fasta" ]; then 
        count_seq=$(cat ${workdir}/${name}.result/${name}.fasta | grep -c '^>');
        printf "ND1 ND2 COX1 COX2 ATP6 COX3 ND3 ND4L ND4 ND5 ND6 ND6 CYTB\n" | \
            tr " " "\n" \
            >${workdir}/${name}.tmp/${name}.assembly/all_genes.txt;
        cut -f2 ${workdir}/${name}.tmp/${name}.assembly/work71.hmmtblout.besthit.sim.filtered | \
            sort -u | \
            fgrep -w -v -f - ${workdir}/${name}.tmp/${name}.assembly/all_genes.txt | \
            sort -u >${workdir}/${name}.tmp/${name}.assembly/missing_genes.txt;
        missing=$(cat ${workdir}/${name}.tmp/${name}.assembly/missing_genes.txt | sort -u | sed -e '/^$/d' | tr "\n" " ");
        if [ ${count_seq} -gt 1 ]; then
            mv ${workdir}/${name}.tmp/${name}.assembly/ ${workdir}/;
            cut -d$'\t' -f1,2 ${workdir}/${name}.assembly/work71.hmmtblout.besthit.sim.filtered | \
                sort -k1,1V -k2,2V | \
                awk -F'\t' '!x[$2]++' | \
                awk -F'\t' 'NF>1 {a[$1] = a[$1]" "$2};END{for(i in a)print i" "a[i]}' \
                >${workdir}/${name}.assembly/quick_mode_fa_genes.txt ;
            if [ -e "${workdir}/${name}.tmp/" ]; then
                mkdir -p ${workdir}/OLDER;
                mv ${workdir}/${name}.tmp/ ${workdir}/${name}.result/ ${workdir}/OLDER;
            fi    
            python3 Software/mitoz/MitoZ.py all2 --run_mode 3 \
                --outprefix ${name} --create_config \
                --fastq1 ${fastq1} --fastq2 ${fastq2} \
                --fastq_read_length 150 --insert_size 450 --fq_size 7 \
                --requiring_taxa 'Chordata' \
                --thread_number 32 \
                --quick_mode_seq_file ${workdir}/${name}.assembly/work71.hmmout.fa \
                --quick_mode_fa_genes_file ${workdir}/${name}.assembly/quick_mode_fa_genes.txt  \
                --quick_mode_score_file ${workdir}/${name}.assembly/work71.hmmtblout.besthit.sim.filtered.high_abundance_10.0X.reformat.sorted  \
                --quick_mode_prior_seq_file ${workdir}/${name}.assembly/work71.hmmtblout.besthit.sim.filtered.fai \
                --missing_PCGs ${missing};
        fi
    fi    
    cd ${maindir};

done    


### python3 Software/mitoz/MitoZ.py all2 --run_mode 3 --outprefix S-rosaceus_SEB-73 --create_config --fastq1 Rockfish/Genomes/sequencing/illumina/output/Sebastes_rosaceus/S-rosaceus_SEB-73/corrected/corr.0 --fastq2 Rockfish/Genomes/sequencing/illumina/output/Sebastes_rosaceus/S-rosaceus_SEB-73/corrected/corr.1 --fastq_read_length 150 --insert_size 450 --fq_size 7 --requiring_taxa 'Chordata' --thread_number 32 --quick_mode_seq_file Rockfish/Genomes/sequencing/illumina/output/Sebastes_rosaceus/S-rosaceus_SEB-73/mitogenome/S-rosaceus_SEB-73.assembly/work71.hmmout.fa --quick_mode_fa_genes_file Rockfish/Genomes/sequencing/illumina/output/Sebastes_rosaceus/S-rosaceus_SEB-73/mitogenome/S-rosaceus_SEB-73.assembly/quick_mode_fa_genes.txt --quick_mode_score_file Rockfish/Genomes/sequencing/illumina/output/Sebastes_rosaceus/S-rosaceus_SEB-73/mitogenome/S-rosaceus_SEB-73.assembly/work71.hmmtblout.besthit.sim.filtered.high_abundance_10.0X.reformat.sorted --quick_mode_prior_seq_file Rockfish/Genomes/sequencing/illumina/output/Sebastes_rosaceus/S-rosaceus_SEB-73/mitogenome/S-rosaceus_SEB-73.assembly/work71.hmmtblout.besthit.sim.filtered.fai --missing_PCGs ND6 && cd Rockfish/Genomes/sequencing/illumina/

