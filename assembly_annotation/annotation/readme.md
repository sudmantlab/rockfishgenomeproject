Annotations were generated for reference and illumina-only genomes

Independent runs, order matters from time to time
1) refgenome_repeatmasker_Snakefile - used to repeatmask reference genomes, same parameters for illumina-only genomes too (tired of copying everything in here)
2) refgenome_starmap_Snakefile - Aliging RNAseq reads to reference genomes, important for annotations 
3) refgenome_funannotate_Snakefile_rnaseq - Protein coding genes annotation using Funannotate (no RNAseq for aleutianus) 
4) refgenome_interproscan_Snakefile - Gene family and functional assignments using coding protein sequences
5) refgenome_noncoding_Snakefile - Noncoding tRNA, mRNA, rRNA, snoRNA annotation using Infernal and tRNAscan-SE
6) refgenome_liftoff_Snakefile_illumina_large -  Lifting annotations from reference genomes to reference based scaffolded illumina-only assemblies 

