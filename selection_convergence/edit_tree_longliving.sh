#!/usr/bin/bash

set -e -o pipefail

file1=$@

if [[ "$#" -ne 1 || ! -f "$file1" ]]; then
    echo "Need an existing input file"
    exit 1;
fi

sed -e 's/Sebastes_aleutianus/Sebastes_aleutianus{Foreground}/' -e 's/Sebastes_alutus/Sebastes_alutus{Foreground}/'  -e 's/Sebastes_aurora/Sebastes_aurora{Foreground}/'  -e 's/Sebastes_babcocki/Sebastes_babcocki{Foreground}/'  -e 's/Sebastes_borealis/Sebastes_borealis{Foreground}/'  -e 's/Sebastes_crameri/Sebastes_crameri{Foreground}/'  -e 's/Sebastes_nigrocinctus/Sebastes_nigrocinctus{Foreground}/'  -e 's/Sebastes_ruberrimus/Sebastes_ruberrimus{Foreground}/'  -e 's/Sebastolobus_alascanus/Sebastolobus_alascanus{Foreground}/' $file1


