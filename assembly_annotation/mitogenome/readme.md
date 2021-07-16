De novo assembly and annotation of mitogenomes from Illumina PE data

Run as follows
1) mitoz_Snakefile - running assembly using MitoZ
2) only if assembly fails - repeat_mitoz_individuall_ifneeded.sh
3) mitos_annotation_part1.sh - annotates
4) mitos_annotation_part2.sh - circularity problem is resolved i.e. mitogenome is ordered by tRNA first and Dloop at the end based on part1 annotation, then reannotated
5) extract_mitogenes.sh - extracting individual genes
