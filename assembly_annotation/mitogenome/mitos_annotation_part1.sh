#!/usr/bin/bash

set -e -o pipefail

# Run on each mitoassembly separately

\ls -d /global/scratch2/Rockfish/Genomes/sequencing/illumina/output/*/*/mitogenome/*result/*.fasta | grep -v 'work' |
    while read fasta1; do
        working=$(echo $fasta1 | sed -e 's/\/mitogenome\/.*/\/mitogenome\//');
        mkdir -p ${working}/MITOS;
        echo "Working at ${working}";
        runmitos.py -i ${fasta1} -o ${working}/MITOS \
            --refdir /global/scratch2/miniconda3/envs/mitos/lib/python2.7/site-packages/mitos/data/refseq81m \
            -c 2 --evalue 1e-6 --best --ncev 1e-3;
    done    


