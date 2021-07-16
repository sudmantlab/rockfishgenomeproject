#!/bin/bash

set -e 

##### Modified for local runs from VGP_assembly github - https://github.com/VGP/vgp-assembly

##### RUN AS - bash run_bionano_solve.sh bionano_config_DLE1.xml DLE1 sample_Saphyr_DLE1_3450142.cmap genome.purged.t1.frby.2.fasta genome_bionano

# TBH we didn't need bionano, our genome was small (<1Gb), and scaffolding with HiC alone was great - well whatever dude/duditt

CONFIG=$1
ENZYME=$2
BMAP=$3
ASM=$4
name=$5
tools=/global/scratch2/Software/bionano_solve/Solve3.4_06042019a

RefAligner=$tools/RefAligner/8949.9232rel/avx/RefAligner

module load gcc gsl openmpi

echo "perl $tools/HybridScaffold/06042019/hybridScaffold.pl \
     -n $ASM -b $BMAP -c $CONFIG -r $RefAligner -B 2 -N 2 -f -o $PWD/$name"

/global/scratch2/Software/perl-5.14.4/bin/perl $tools/HybridScaffold/06042019/hybridScaffold.pl \
    -n $ASM -b $BMAP -c $CONFIG -r $RefAligner -B 2 -N 2 -f -o $PWD/$name

