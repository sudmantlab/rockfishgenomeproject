Scripts used for associating lifespan and substitution rate.

1) Snakefile_tree_dnds - Creating dNdS trees using Hyphy (Need alignment, species tree, Hyphy2 installed). Activate hyphy2.yaml for usage
2) Snakefile_tree_dnds-pairdist - Calculating pairwise substitution rates. Activate R.yaml for usage (R.yaml is heavy and has multiple other packages used in other workflows too beyond this particular analysis)
3) substitution_rates.Rmd - R notebook for code to compare substitutions rates of coding sequences based on tree (uses nuclear_speciespairs_dnds_calculated.txt file which is the collection of dN and dS values between species pairs)
The dN, dS trees were produced from Hyphy package using their script dNdS.bf
The Poisson regression tests has been taken from Hua et al. 2015
