Annotations were generated for reference and illumina-only genomes

Independent runs, order matters from time to time
1) refgenome_repeatmasker_Snakefile - used to repeatmask reference genomes, same parameters for illumina-only genomes too (tired of copying everything in here)
2) refgenome_starmap_Snakefile - Aliging RNAseq reads to reference genomes, important for annotations 
3) refgenome_funannotate_Snakefile_rnaseq - Protein coding genes annotation using Funannotate (no RNAseq for aleutianus) 
4) refgenome_interproscan_Snakefile - Gene family and functional assignments using coding protein sequences
5) refgenome_noncoding_Snakefile - Noncoding tRNA, mRNA, rRNA, snoRNA annotation using Infernal and tRNAscan-SE
6) refgenome_liftoff_Snakefile_illumina_large -  Lifting annotations from reference genomes to reference based scaffolded illumina-only assemblies 
7) probable.orthologs.tsv.zip - Orthologs between the 88 species. The reference genomes were used to identify orthologs using Protheinortho6 tool (refer to Ortholog section in the manuscript supplement), and the liftoff annotations were used for the Illumina-only genomes. Post-processing was performed based on alignment quality and divergence (refer to selection_convergence directory)
8) The annotation files can be found on Zenodo - https://zenodo.org/record/5534983
