#!/usr/bin/env Rscript

if (!require("dplyr")) install.packages("dplyr")
if (!require("ape")) install.packages("ape")
if (!require("assertthat")) install.packages("assertthat")

library(ape)
library(assertthat)
library(RERconverge)
library(data.table)
library(rlist)
library(stringr)
library(reshape2)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
    stop("Input arguments are species_tree weights_file work_directory : Now given only ",
    paste0(length(args)), call.=FALSE)
}

seb_tree_77 <- ape::read.tree(args[1], keep.multi = F)
seb_tree_77_path <- args[1]
given_tree_weight <- args[2]
outdirectory = args[3]

weight_value <- scan(given_tree_weight, character(), quote = "")
weight_value <- weight_value %>%
    unlist %>% unname %>% as.vector

setwd(outdirectory)
print(outdirectory)

pruned_ages <- fread(input = "rockfish_ages.txt.filt")
row.names(pruned_ages) <- pruned_ages$species

pruned_ages <- pruned_ages %>%
  dplyr::select(species, lifespan) %>%
  mutate_each(funs(as.numeric), lifespan)

pruned_ages <- pruned_ages[match(seb_tree_77$tip.label, pruned_ages$species),]

age_list <- pruned_ages$lifespan
names(age_list) <- pruned_ages$species

age_list_1 <- age_list[!is.na(age_list)]

estimatePhangornTreeAll(alndir="./select_filter_trim80/",
                        pattern=".FAA", treefile="./tree.nwk.2", output.file="./gene_trees.tree", 
                        type="AA", format="fasta"); 



