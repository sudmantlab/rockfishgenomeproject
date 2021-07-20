#!/use/bin/env Rscript

if (!require("dplyr")) install.packages("dplyr")
if (!require("ape")) install.packages("ape")
if (!require("assertthat")) install.packages("assertthat")

library(ape)
library(dplyr)
library(assertthat)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Input arguments are species_tree fasta_indexfile output_treefile : Now given only ", 
       paste0(length(args)), call.=FALSE)
}

species_tree <- ape::read.tree(args[1], keep.multi = F)
species_index <- read.csv2(args[2], header = F, sep = "\t",
                           blank.lines.skip = TRUE)
outtree = args[3]

species_names <- species_index %>%
  dplyr::select(V1) %>%
  as.vector() %>% unlist() %>% unname() %>% levels()

all_species <- species_tree$tip.label
to_remove <- setdiff(all_species, species_names)

out_tree <- ape::drop.tip(species_tree, to_remove)
ape::write.tree(out_tree, file=outtree)

