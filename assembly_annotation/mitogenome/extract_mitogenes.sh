#!/usr/bin/bash

set -e

module load seqtk samtools

    grep -v -e O[HL] -e rrn -e trn mitogenome/MITOS/result.bed |
    sed -e '/^$/d' | 
    awk -F'\t' '{ print $4}' >mitogenome/list_mitogenes.txt

\ls -d /global/scratch2/Rockfish/Genomes/sequencing/illumina/output/*/*/mitogenome/*.result/*.fasta | grep -v 'work' |
    sed -e 's/\/mitogenome\/.*//' -e 's/.*\/output\///' -e 's/\//\t/' |
    sort -k1,1 -k2,2 | awk -F'\t' '!x[$1]++ {print $1"-"NR"\t"$2}' >list_names.txt

if [ -e "list_problematic.txt" ]; then 
    rm list_problematic.txt; 
fi

count=1
\ls -d /global/scratch2/Rockfish/Genomes/sequencing/illumina/output/*/*/mitogenome/*.result/*.fasta | 
    grep -v 'work' |
        while read fasta1; do
            working=$(echo $fasta1 | sed -e 's/\/mitogenome\/.*/\/mitogenome\//');
            speciesname=$(echo $fasta1 | sed -e 's/\/mitogenome\/.*//' -e 's/.*\/output\///' -e 's/\/.*//');
            samplename=$(echo $fasta1 | sed -e 's/\/mitogenome\/.*//' -e 's/.*\/output\///' -e 's/.*\///');

            if [ ! -e "${working}/MITOS/CDS/ALL.FAS" ]; then
	            samtools faidx $fasta1;
	            size=$(awk '{print $2}' ${fasta1}.fai);

        	    printf "\tProcessing - ${speciesname}\n";
	            mkdir -p ${working}/MITOS/CDS;
        	    touch ${working}/MITOS/CDS/${speciesname}.fas

	            cat list_mitogenes.txt | while read genename; do
        	        printf "\tExtracting - ${genename}\n";
		        grep -w "${genename}" ${working}/MITOS/result.bed >${working}/MITOS/CDS/${genename}.bed.edit
        		awk -v size="$size" '{if($3<$2) {print $1"\t"$2"\t"size"\t"$1"\t"0"\t"$3} else print $1"\t"$2"\t"$3}' \
                	    ${working}/MITOS/CDS/${genename}.bed.edit >${working}/MITOS/CDS/${genename}.bed        ### Circularity Problem
               		rm ${working}/MITOS/CDS/${genename}.bed.edit
	                seqtk subseq ${fasta1} ${working}/MITOS/CDS/${genename}.bed |
        	            grep -v '^>' |
                	    sed -e "1i >${speciesname} ${genename}" >${working}/MITOS/CDS/${genename}.FAS
                	cat ${working}/MITOS/CDS/${genename}.FAS | grep -v '^>' >>${working}/MITOS/CDS/${speciesname}.fas    
        	    done

            	printf ">${speciesname}-${count}\n" >${working}/MITOS/CDS/ALL.FAS
	        cat ${working}/MITOS/CDS/${speciesname}.fas | 
        	        tr "\n" " " | 
                	sed -e 's/ //g' -e 's/>/\n>/' -e '/^$/d'  >>${working}/MITOS/CDS/ALL.FAS
                rm ${working}/MITOS/CDS/${speciesname}.fas
            fi    
            
	    problem_count=$(ls -l ${working}/MITOS/CDS/[a-z]*.FAS | awk '$5<50' | wc -l)    

            if [ "${problem_count}" -gt 0 ]; then	# Problematic assemblies
                printf "${speciesname}\t${fasta1}" >>list_problematic.txt
            fi

            count=$((count+1))
done

if [ -e "list_problematic.txt" ]; then 
    problems=$(wc -l list_problematic.txt);
    print "\n\nPLEASE CHECK FOR ${problems} ERROR SAMPLES :  - list_problematic.txt\n";
fi



