id="sebasto_sebaste_acti"
rm ${id}_per_gene_aligned_filt_ts_distfilt.treefile; 
for i in `ls ${id}_per_gene_aligned_filt_ts_distfilt_fasta | sed s/_NT.fasta//g`; do cat ${id}_per_gene_trees_unrooted/$i.treefile >> ${id}_per_gene_aligned_filt_ts_distfilt.treefile; done
nw_ed ${id}_per_gene_aligned_filt_ts_distfilt.treefile  'i & b<=10' o >  ${id}_per_gene_aligned_filt_ts_distfilt_bs10.treefile
java -jar /home/owens/bin/ASTRAL-5.7.1/astral.5.7.1.jar -i ${id}_per_gene_aligned_filt_ts_distfilt_bs10.treefile -o ${id}_per_gene_aligned_filt_ts_distfilt_bs10.astral.tr -t 3
bash ../../bin/rename_reverse.sh ${id}_per_gene_aligned_filt_ts_distfilt_bs10.astral.tr meta/sample_renames.txt > ${id}_per_gene_aligned_filt_ts_distfilt_bs10.astral.renamed.tr
