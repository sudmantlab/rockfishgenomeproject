#!/usr/bin/bash

set -e -o pipefail

file1=$@

if [[ "$#" -ne 1 || ! -f "$file1" ]]; then
    echo "Need an existing input file"
    exit 1;
fi    

sed  -e 's/Sebastes_dalli/Sebastes_dalli{Foreground}/'  -e 's/Sebastes_hopkinsi/Sebastes_hopkinsi{Foreground}/'  -e 's/Sebastes_inermis/Sebastes_inermis{Foreground}/'  -e 's/Sebastes_minor/Sebastes_minor{Foreground}/'  -e 's/Sebastes_semicinctus/Sebastes_semicinctus{Foreground}/'  -e 's/Sebastes_steindachneri/Sebastes_steindachneri{Foreground}/'  $file1

