#!/usr/bin/bash

set -e

conda activate pb-assembly
mkdir -p run1 && cd run1; 

# We are using these config files are optimised to run locally on a high memory (1.5TB RAM) machine with 32 threads - these are not for compute farms 
# YESSSSS - I was able to tame the falcon beast locally after many works of burning our allocation units - dont try this at home/work/office, just use my configs

fc_run.py fc_run_aleutianus.cfg 1>>LOG 2>>LOG ;
fc_unzip.py fc_unzip_aleutianus.cfg 1>LOG2 2>>LOG2

