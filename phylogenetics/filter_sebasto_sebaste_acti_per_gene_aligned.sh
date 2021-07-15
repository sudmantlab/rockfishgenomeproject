rm sebasto_sebaste_acti_per_gene_aligned_filt_bases.txt;
for i in `ls sebasto_sebaste_acti_per_gene_aligned | grep NT | sed s/.FNA//g`; do  /home/owens/bin/phyx/src/pxclsq -s sebasto_sebaste_acti_per_gene_aligned/$i.FNA -p 0.3 | sed s/\!/\-/g > sebasto_sebaste_acti_per_gene_aligned_filt/$i.fas; bases=$(tail -n 1 sebasto_sebaste_acti_per_gene_aligned_filt/$i.fas | wc -m); echo "$i $bases" >> sebasto_sebaste_acti_per_gene_aligned_filt_bases.txt; done
for i in `awk '($2 < 1000)' sebasto_sebaste_acti_per_gene_aligned_filt_bases.txt | cut -f 1 -d " "`; do  rm sebasto_sebaste_acti_per_gene_aligned_filt/$i.fas; done