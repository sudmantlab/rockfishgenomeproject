Scripts to run RERconverge analysis and creating the ternary plots

1) astral_samples_trees.nwk - all the tree topologies from ASTRAL 
2) Snakefile_rerconverge - runs RERconverge on a directory with subdirectories of different tree topologies, the trait used was max lifespan (rockfish_ages.txt.filt). Note that some topologies might need to be discarded if convergence cannot be achieved for particular gene trees in the estimatePhangornTree step of RERconverge

RERconverge can be repeated for other traits using LH_rockfish_34.txt

3) plot_RER_converge_analyses_final - directory with Rnotebook to run Ternary plot analyses
