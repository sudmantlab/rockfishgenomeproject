for n in `seq 711`; do file=$(ls sebasto_sebaste_acti_per_gene_aligned_filt | head -n $n | tail -n 1 | sed s/.fas//g); mkdir -p sebasto_sebaste_acti_per_gene_aligned_filt_ts/$file; ln -s -f /home/owens/working/Sebastes/data/trees/sebasto_sebaste_acti_per_gene_aligned_filt/$file.fas sebasto_sebaste_acti_per_gene_aligned_filt_ts/$file/input.fasta; cat sebasto_sebaste_acti_per_gene_aligned_filt_genetrees/loci.treefile | head -n $n | tail -n 1 > sebasto_sebaste_acti_per_gene_aligned_filt_ts/$file/input.tree; done

run_treeshrink.py -i sebasto_sebaste_acti_per_gene_aligned_filt_ts/ -t input.tree -a input.fasta 1> sebasto_sebaste_acti_per_gene_aligned_filt_ts/log.txt

for gene in `ls sebasto_sebaste_acti_per_gene_aligned_filt_ts | grep NT`; do ln -s /home/owens/working/Sebastes/data/trees/sebasto_sebaste_acti_per_gene_aligned_filt_ts/$gene/output.fasta sebasto_sebaste_acti_per_gene_aligned_filt_ts_fasta/$gene.fasta; done


