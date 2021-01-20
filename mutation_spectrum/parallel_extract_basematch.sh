chr=$1

cat /global/scratch2/rohitkolora/Rockfish/Data/Greg/snp_calling/aleutianus_scaffolded/SPECIES/align.Sebastes_aleutianus.SPECIES.maf | grep -B 1 -A 1 "^s $chr" | perl /global/scratch/gregoryowens/sebastes/bin/extracted_paired_maf_positions.pl /global/scratch2/rohitkolora/Rockfish/Data/Greg/snp_calling/aleutianus_scaffolded/SPECIES/Illumina/map.SPECIES.SAMPLE.hetsites.filtered.txt SPECIES > /global/scratch2/rohitkolora/Rockfish/Data/Greg/snp_calling/all_aleutianus_hetsites/SPECIES/map.$chr.SPECIES.SAMPLE.tmpset.sites.txt
