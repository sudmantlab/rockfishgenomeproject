We de novo assemble the illumina PE data

Run as follows -
1) masurca_Snakefile - for whole genome de novo assembly using MASURCA
The above also annotates the genome with BRAKER2, but we dont use these anyway - have fun reading the code
2) busco_masurca_Snakefile - for BUSCO score check
3) refgenome_ragtag_Snakefile - for annotating coding genes based on reference genomes 
