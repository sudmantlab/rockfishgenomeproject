#!/usr/bin/bash

#set -e -o pipefail

outname=$4
outname_bf=${outname}.${BASHPID}.bf

if [ $# -ne 4 ]; then
    echo "  Need the genetic_code alignment_file.nex tree_file.nwk"
    printf "\tbash ~/RGP/scripts/selection/hyphy_dnds.run.sh Universal infile.nex infile.nwk outname \n\n"
    exit 1
fi

bash ~/RGP/scripts/selection/hyphy_dnds.inputcreate.sh $1 $2 $3 >${outname_bf}

hyphy ${outname_bf} >${outname}.LOG

grep -m 1 -A1 -e '^dN tree:' ${outname}.LOG | sed -e 's/ //g' -e 's/\t//g' -e 's/$/;/' | awk 'NR==2' >${outname}.dN.nwk
grep -m 1 -A1 -e '^dS tree:' ${outname}.LOG | sed -e 's/ //g' -e 's/\t//g' -e 's/$/;/' | awk 'NR==2' >${outname}.dS.nwk
grep -m 1 'Total nucleotide sites :' ${outname}.LOG | sed -e 's/Total nucleotide sites ://' >${outname}.length
grep -m 1 'Synonymous  sites ' ${outname}.LOG | sed -e 's/.*://' >${outname}.synsites
grep -m 1 'Nonsynonymous  sites ' ${outname}.LOG | sed -e 's/.*://' >${outname}.nonsynsites
