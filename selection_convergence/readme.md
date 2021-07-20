Analyses for selection (Codeml & Hyphy - absrel, busted etc) and convergence

1) Snakefile_hyphy_aln - creating alignments of orthogroups for downstream analyses (each subdirectory with sample name has individual fasta and related files)
2) Snakefile_hyphy_tree - creating tree files by retaining needed tips/species
3) misc_hyphy_commands.txt - commands used to run different tests for selection using HyPhy
4) Snakefile_codeml - running codeml for long lived species
5) codeml_calculate_pval.sh - using codeml output for performing a Chi-square test to test for signifance in selection (useful to get important sites)
6) Snakefile_converged_aa - looking at AA convergence in protein sequence alignments for long/short lived species
