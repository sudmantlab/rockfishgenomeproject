for n in `seq 659`; do file=$(ls all_species_acti_per_gene_aligned_filt | head -n $n | tail -n 1 | sed s/.fas//g); mkdir -p all_species_acti_per_gene_aligned_filt_ts/$file; ln -s -f /home/owens/working/Sebastes/data/trees/all_species_acti_per_gene_aligned_filt/$file.fas all_species_acti_per_gene_aligned_filt_ts/$file/input.fasta; cat all_species_acti_per_gene_aligned_filt_genetrees/loci.treefile | head -n $n | tail -n 1 > all_species_acti_per_gene_aligned_filt_ts/$file/input.tree; done

run_treeshrink.py -i all_species_acti_per_gene_aligned_filt_ts/ -t input.tree -a input.fasta 1> all_species_acti_per_gene_aligned_filt_ts/log.txt

for gene in `ls all_species_acti_per_gene_aligned_filt_ts | grep NT`; do ln -s /home/owens/working/Sebastes/data/trees/all_species_acti_per_gene_aligned_filt_ts/$gene/output.fasta all_species_acti_per_gene_aligned_filt_ts_fasta/$gene.fasta; done


