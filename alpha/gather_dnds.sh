base_dir="/global/scratch2/rohitkolora/Rockfish/Genomes/orthologs/Lifted/Sebastes_55/select_sequence_files"

rm seq_names.txt;
rm dnds_values.txt;
for orthogroup in `ls $base_dir | grep OrthoGroup `;
do
cat $base_dir/$orthogroup/seqs.names | sed s/^/$orthogroup\\t/g >> seq_names.txt
#echo >> seq_names.txt

cat $base_dir/$orthogroup/Filter/dNdS_orf/$orthogroup.dS.nwk | tr ',' '\n' | tr '(' '\n' | tr ')' '\n' | grep Seb | sed s/\:/\\t/g | sed s/^/$orthogroup\\tdS\\t/g >> dnds_values.txt
 
#echo >> dnds_values.txt

cat $base_dir/$orthogroup/Filter/dNdS_orf/$orthogroup.dN.nwk | tr ',' '\n' | tr '(' '\n' | tr ')' '\n' | grep Seb | sed s/\:/\\t/g | sed s/^/$orthogroup\\tdN\\t/g >> dnds_values.txt

#echo >> dnds_values.txt

done
