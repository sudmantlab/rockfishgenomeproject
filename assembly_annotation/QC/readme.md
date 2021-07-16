Quality checks were performed for both reference assemblies and illumina only assemblies

Check as following
1) refgenome_busconew_Snakefile && refgenome_busco4peps_Snakefile - BUSCO for reference genomes with different modes
2) busco_masurca_Snakefile - BUSCO for illumina only genomes (BUSCO and its alterations, use diff modes, diff databases, diff versions, diff blah blah)
3) refgenome_meryl_Snakefile - creating meryl databases for the reference genome assemblies
4) refgenome_merqury_Snakefile - quality check using Merqury tool (thanks Arang!!!)

