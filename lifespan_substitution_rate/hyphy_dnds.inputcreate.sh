#!/usr/bin/bash

HYPHY_LIB_DIRECTORY="/global/scratch2/miniconda3/envs/hyphy/lib/hyphy/"

genetic_code=$1
alignment_file=$2
tree_file=$3

cat << EOF
inputRedirect = {};
inputRedirect["01"]="${genetic_code}"; // genetic code    -  Universal Vertebrate-mtDNA
inputRedirect["02"]="Custom"; // for input without example
inputRedirect["03"]="${alignment_file}"; // codon data      -  nexus alignment file
inputRedirect["04"]="${tree_file}"; // tree            -   newick tree file
inputRedirect["05"]=""; 
//inputRedirect["05"]="1"; // Test for selection on a branch    -   1=ALL branches to test
//inputRedirect["06"]=""; // complete selection

ExecuteAFile ("${HYPHY_LIB_DIRECTORY}"+"TemplateBatchFiles/Miscellaneous/phylohandbook/dSdN.bf", inputRedirect);
EOF


