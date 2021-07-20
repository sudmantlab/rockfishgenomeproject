## Scripts for calculating the mutation spectrum and it's association with lifespan.

Description:

- compare_hetsites_to_all_others.pl: Checks heterozygous sites to see if one allele is unique to that particular species and the site is genotyped in sufficient other species.
- extracted_paired_maf_positions.pl: Pulls out matching position and base call between S. aleutianus and another species genome.
- extracted_paired_maf_positions_otherbase.pl: Extracts other positions used for checking allelic state in other species.
- filter_mappability_bed.pl: From output of genmap, creates a bedfile of maximum mappability regions.
- filter_site_spacing.pl: Filters heterozygous sites that are beside other heterozygous sites.
- gather_hetsites.txt: Gathers together big lists of heterozygous sites used so we have all possible sites in a single file.
- get_aleutianus_bases.pl: Grabs surrounding bases of heterozygous positions for S. aleutianus
- get_surrounding_bases.pl: Grabs surrounding bases of heterozygous positions for other species.
- mutation_spectra_scaf.R: Exploration of mutation spectrum
- parallel_extract_basematch.sh: Parallel running of extracted_paired_maf_positions.pl
- parallel_extract_basematch_all.sh: Parallel running of extracted_paired_maf_positions_otherbase.pl
- scafgenome_align_pt1_Snakefile: Aligning all species genomes
- scafgenome_align_pt2_Snakefile: Realigning, filtering, calling variants and measuring mappability
- scafgenome_align_pt3_Snakefile: Measuring heterozygousity.
