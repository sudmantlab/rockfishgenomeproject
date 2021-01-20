chr=$1

cat /global/scratch2/rohitkolora/Rockfish/Data/Greg/snp_calling/aleutianus_scaffolded/SPECIES/align.Sebastes_aleutianus.SPECIES.maf | grep -B 1 -A 1 "^s $chr" | perl /global/scratch/gregoryowens/sebastes/bin/extracted_paired_maf_positions_otherbase.pl /global/scratch2/rohitkolora/Rockfish/Data/Greg/snp_calling/all_aleutianus_hetsites/ALL_sites.txt.gz $chr > /global/scratch2/rohitkolora/Rockfish/Data/Greg/snp_calling/all_aleutianus_hetsites/SPECIES/map.$chr.SPECIES.tmpset.sites.txt
