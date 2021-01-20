rm all_species_acti_per_gene_aligned_filt_bases.txt;
for i in `ls all_species_acti_per_gene_aligned | grep _NT | sed s/.FNA//g`; do  /home/owens/bin/phyx/src/pxclsq -s all_species_acti_per_gene_aligned/$i.FNA -p 0.3 | sed s/\!/\-/g > all_species_acti_per_gene_aligned_filt/$i.fas; bases=$(tail -n 1 all_species_acti_per_gene_aligned_filt/$i.fas | wc -m); echo "$i $bases" >> all_species_acti_per_gene_aligned_filt_bases.txt; done
for i in `awk '($2 < 1000)' all_species_acti_per_gene_aligned_filt_bases.txt | cut -f 1 -d " "`; do  rm all_species_acti_per_gene_aligned_filt/$i.fas; done
