## Scripts for calculating effective population size and it's association with lifespan.

Description:
- filter_mappability_bed.pl: Takes a file of mappability scores from genmap and creates a bed file of regions with maximum mappability
- msmc_plots.R: Plotting msmc results
- msmc_vs_lifespan.R: Comparison of recent Ne with maximum lifespan
- prep_msmc_Snakefile: Mask and prep vcf for use in msmc
- refgenome_illuminamap_Snakefile: Mapping and variant calling for individual samples
- refgenome_neprep_perscaffold_Snakefile: Mapping and prepping for msmc, but run per scaffold to improve speed
- run_msmc_scafgenome_Snakefile: Run msmc, excluding shorter scaffolds
