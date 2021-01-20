base_dir="/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/select_sequence_files"

rm dn_lengths.txt;
rm ds_lengths.txt;
for orthogroup in `ls $base_dir | grep OrthoGroup `;
do
if [ ! -f "$base_dir/$orthogroup/Filter/dNdS_orf/$orthogroup.nonsynsites" ]
then
continue;
fi;
cat $base_dir/$orthogroup/Filter/dNdS_orf/$orthogroup.nonsynsites| sed s/^/$orthogroup\\t/g >> dn_lengths.txt

cat $base_dir/$orthogroup/Filter/dNdS_orf/$orthogroup.synsites| sed s/^/$orthogroup\\t/g >> ds_lengths.txt

done
