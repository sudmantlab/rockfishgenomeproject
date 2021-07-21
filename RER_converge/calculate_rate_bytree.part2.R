#!/usr/bin/env Rscript

if (!require("dplyr")) install.packages("dplyr")
if (!require("ape")) install.packages("ape")
if (!require("assertthat")) install.packages("assertthat")

options(show.error.locations = TRUE)

library(ape)
library(assertthat)
library(RERconverge)
library(data.table)
library(rlist)
library(stringr)
library(reshape2)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=4) {
    stop("Input arguments are species_tree weights_file work_directory : Now given only ",
    paste0(length(args)), call.=FALSE)
}

seb_tree_77 <- ape::read.tree(args[1], keep.multi = F)
given_tree_weight <- args[2]
outdirectory = args[3]
sample_number = args[4]

weight_value <- scan(given_tree_weight, character(), quote = "")
weight_value <- weight_value %>%
    unlist %>% unname %>% as.vector %>% head(n=1L)

setwd(outdirectory)

pruned_ages <- fread(input = "rockfish_ages.txt.filt")
row.names(pruned_ages) <- pruned_ages$species

pruned_ages <- pruned_ages %>%
  dplyr::select(species, lifespan) %>%
  mutate_each(funs(as.numeric), lifespan)

pruned_ages <- pruned_ages[match(seb_tree_77$tip.label, pruned_ages$species),]

age_list <- pruned_ages$lifespan
names(age_list) <- pruned_ages$species

age_list_1 <- age_list[!is.na(age_list)]

select_all_trees <- RERconverge::readTrees("gene_trees.tree")

age_charpaths=char2Paths(age_list_1, select_all_trees)

pdf("longevity.plot_residuals.pdf")
residuals_all_genes <- getAllResiduals(select_all_trees, 
                                       useSpecies=names(age_list_1), min.sp = 15,
                                       transform = "sqrt", weighted = T, scale = T, plot = T)
dev.off()

cor_age_all_genes=correlateWithContinuousPhenotype(residuals_all_genes,
                                         age_charpaths, min.sp=15,
                                         winsorizeRER = 3, winsorizetrait = 3)
cor_age_all_genes$Genes <- rownames(cor_age_all_genes)
cor_age_all_genes$Treeweight <- c(weight_value)
cor_age_all_genes$Sample <- c(sample_number)

write.table(cor_age_all_genes,
            file="longevity.correlation_ages.txt",
            quote = F, sep = "\t", row.names = F, col.names = T)

