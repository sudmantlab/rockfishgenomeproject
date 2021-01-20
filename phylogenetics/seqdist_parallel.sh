#i=$1
#if [ ! -f sebasto_sebaste_acti_per_gene_aligned_filt_ts_seqdist/$i.seqdist ]
#then
#/home/owens/working/Sebastes/bin/distats/distats.pl --print_dist_matrix sebasto_sebaste_acti_per_gene_aligned_filt_ts_fasta/$i.fasta ${i}_dist;
#mv distance_matrix_${i}_dist sebasto_sebaste_acti_per_gene_aligned_filt_ts_seqdist/$i.seqdist; rm ${i}_dist;
#fi

i=$1
if [ ! -f all_species_acti_per_gene_aligned_filt_ts_seqdist/$i.seqdist ]
then
/home/owens/working/Sebastes/bin/distats/distats.pl --print_dist_matrix all_species_acti_per_gene_aligned_filt_ts_fasta/$i.fasta ${i}_dist;
mv distance_matrix_${i}_dist all_species_acti_per_gene_aligned_filt_ts_seqdist/$i.seqdist; rm ${i}_dist;
fi
