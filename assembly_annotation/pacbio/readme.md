For de novo whole genome assembly and subsequent processing, here we currently use PacBio sequel-I, Illumina, HiC and Bionano data

Run in this order
# Polishing snd bionano steps were added & modified based on the VGP_assembly protocol - https://github.com/VGP/vgp-assembly
1) extract_pacbiorun_tar.sh - extracting sequencing data from sequencing center raw files 
2) extract_subreads_bam2fastq.sh - extracting data 
3) misc_merge_fastq.sh - combining things
4) run_falcon.sh - produces primary and alternate assemblies
5) run_purging.sh - for primary assembly
6) run_polish_arrow.sh - for both primary and alternate
7) run_polish_freebayes.sh - for both primary and alternate
8) run_bionano.sh - if present
9) run_scaffolding_hic_salsa_dnase.sh || run_scaffolding_hic_juicer-3ddna.sh
10) After manual review using juicebox for the scaffolds convert using this command : py2.7 juicer/1.6.2/misc/juicebox_assembly_converter.py -a ./EDIT_GENOME/scaffolds_FINAL.0.review.assembly -f salsa/scaffolds/scaffolds_FINAL.fasta -v -p -c ./EDIT_GENOME/genome_falcon_convert 
11) rockfish_named_chromosomes.stranded.txt - The chromosome numbering and the strand direction wrt S. aleutianus used for detecting structural variations across rockfish

