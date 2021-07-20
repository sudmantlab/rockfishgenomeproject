#!/usr/bin/bash

set -e -o pipefail

file1=$@

if [[ "$#" -ne 1 || ! -f "$file1" ]]; then
    echo "Need an existing input file"
    exit 1;
fi

sed -e 's/Sebastes_aleutianus/Sebastes_aleutianus #1/' -e 's/Sebastes_alutus/Sebastes_alutus #1/'  -e 's/Sebastes_aurora/Sebastes_aurora #1/'  -e 's/Sebastes_babcocki/Sebastes_babcocki #1/'  -e 's/Sebastes_borealis/Sebastes_borealis #1/'  -e 's/Sebastes_crameri/Sebastes_crameri #1/'  -e 's/Sebastes_nigrocinctus/Sebastes_nigrocinctus #1/'  -e 's/Sebastes_ruberrimus/Sebastes_ruberrimus #1/'  -e 's/Sebastolobus_alascanus/Sebastolobus_alascanus #1/' $file1


