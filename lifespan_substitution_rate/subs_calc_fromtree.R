#!/use/bin/env Rscript

if (!require("dplyr")) install.packages("dplyr")
if (!require("ape")) install.packages("ape")
if (!require("assertthat")) install.packages("assertthat")

library(ape)
library(reshape2)
library(dplyr)
library(assertthat)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=5) {
  stop("Input arguments are -> species_tree num_sites quoted_mut_type quoted_name outfile_name : Now given only ", 
       paste0(length(args)), call.=FALSE)
}

dist_tips_allpairs <- function(input_tree, numsites, quoted_mut_type, quoted_name) {
  group_name <- as.character(quoted_name)
  mut_type <- as.character(quoted_mut_type)
  num_total_sites <- as.numeric(numsites)
  input_tree_mrca <- setNames(melt(mrca(input_tree)), c('species1', 'species2', 'node'))  #mrca node finder pairs
  input_tree_tips <- as.data.frame(input_tree$tip.label)         #species tip labels
  colnames(input_tree_tips) <- c("species")
  input_tree_tips$tipnum <- rownames(input_tree_tips)
  input_tree_expand <- input_tree_mrca %>%
    left_join(input_tree_tips, by = c("species1" = "species")) %>%
    left_join(input_tree_tips, by = c("species2" = "species"))         #merging the tips with mrca nodes
  input <- as.data.frame(input_tree_expand)
  output_df_pair_anc <- data.frame(species1=character(), species2=character(),
                              tipnum.sp1=character(), tipnum.sp2=character(), 
                              mut_type=character(), totsites=double(),
                              node=character(), name=character(),
                              subrate.sp1=double(), subrate.sp2=double(),
                              numsubs.sp1=double(), numsubs.sp2=double(),
                              stringsAsFactors=FALSE)       #dataframe for species pairs, tipnums, nodenums, subrates
  for (row_num in 1:nrow(input)) {
    species1_name <- input[row_num, "species1"]         #species names
    species2_name  <- input[row_num, "species2"]
    species1_tip <- input[row_num, "tipnum.x"]         #species tip numbers
    species2_tip <- input[row_num, "tipnum.y"]
    mrca_node <- input[row_num, "node"]
    if ((species1_name != species2_name)) {        # removing same species
      dist_species1 <- as.numeric( dist.nodes(input_tree)[species1_tip, mrca_node] )  #sub.rates for species to node
      dist_species2 <- as.numeric( dist.nodes(input_tree)[species2_tip, mrca_node] ) 
      sites_species1 <- dist_species1 * num_total_sites 
      sites_species1 <- round(sites_species1, 0)
      sites_species2 <-  dist_species2 * num_total_sites
      sites_species2 <- round(sites_species2, 0)
      row = data.frame(name = group_name, mut_type = mut_type, 
                       totsites = num_total_sites,
                       node = mrca_node,
                       species1 = species1_name, species2 = species2_name,
                       tipnum.sp1 = species1_tip, tipnum.sp2 = species2_tip, 
                       subrate.sp1 = dist_species1, subrate.sp2 = dist_species2,
                       numsubs.sp1 = sites_species1, numsubs.sp2 = sites_species2)
      output_df_pair_anc <- rbind(output_df_pair_anc, row)         
    }
  }
  #output_df_pair_anc
  output_df_pairs <- output_df_pair_anc %>% dplyr::select(-c(node, tipnum.sp1, tipnum.sp2))
  colnames(output_df_pairs) <- c("name", "mut_type", 
                                paste0(mut_type,".sites",sep=""), 
                                "sp1", "sp2", 
                                paste0(mut_type,".subrate.sp1",sep=""),
                                paste0(mut_type,".subrate.sp2",sep=""),
                                paste0(mut_type,".numsubs.sp1",sep=""), 
                                paste0(mut_type,".numsubs.sp2",sep=""))
  output_df_pairs %>% dplyr::select(-c(mut_type))
}

### Reading and outputting
input_tree <- ape::read.tree(args[1], keep.multi = F)
num_sites <- scan(args[2], character(), quote = "")
mut_type = args[3]
mut_type <- as.character(mut_type)
name_group = args[4]
name_group <- as.character(name_group)
outfile_name <- args[5]

output_pairs_subs <- dist_tips_allpairs(input_tree, num_sites, mut_type, name_group)

write.table(output_pairs_subs, outfile_name, 
          row.names=FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

