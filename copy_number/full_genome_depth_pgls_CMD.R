
library(tidyverse)
library(cowplot)
library(nlme)
library(ape)
args = commandArgs(trailingOnly=TRUE)

scaf <- args[1]
#scaf <- "PGA_scaffold_21__11"
sample_depth <- read_tsv("../../data/gene_copies/average_sample_depth.txt")

max_life = read_tsv("../../data/meta/max_life_jul5_2020.txt") %>%
  filter(grepl("Sebastes_",species)) %>%
  filter(!is.na(Max_life))

#updated to new tree March 2021.
tree <- read.tree("../../data/trees/sebasto_sebaste_acti_per_gene_aligned_filt_ts_distfilt_astral/sebasto_sebaste.astral.iqdist.rooted.treefile")




files <- list.files("../../data/gene_copies/aleutianus_genome/full_genome_scan", scaf)

region_depth <- tibble(scaffold=character(),start=numeric(),end=numeric(),
                       gene=character(),total_depth=numeric(),
                       species=character(),sample=character())
for (i in files){
  species <- str_split(i,"\\.")[[1]][2]
  sample <- str_split(i,"\\.")[[1]][3]
  data <- read_tsv(paste0("../../data/gene_copies/aleutianus_genome/full_genome_scan/",i),
                   col_names = c("scaffold","start","end","total_depth"))
  data$sample <- sample
  data$species <- species
  region_depth <- rbind(region_depth,data)
}

region_depth %>%
  inner_join(.,sample_depth) %>%
  filter(species  != "Sebastes_aleutianus") %>%
  mutate(scaled_depth = total_depth/median_depth/100) %>%
  group_by(scaffold,start,end,species) %>%
  summarize(mean_scaled_depth = mean(scaled_depth)) %>%
  ungroup() %>%
  dplyr::select(scaffold, start, end,species, mean_scaled_depth) %>%
  inner_join(max_life) %>%
  filter(species %in% tree$tip.label) -> region_depth_scaled

species_used <- unique(region_depth_scaled$species)
subset_tree <- keep.tip(tree, species_used)

output <- matrix(ncol=4, nrow=length(unique(region_depth_scaled$start)))
count<- 0
for (i in unique(region_depth_scaled$start)){
#for (i in head(unique(region_depth_scaled$start),100)){
  
  count<- count+1
  print(count)
  region_depth_scaled %>%
    filter(start == i)  %>%
    as.data.frame -> test
  if (var(test$mean_scaled_depth) == 0){next;}
  
  row.names(test) <- test$species
  pglsModel <- gls(mean_scaled_depth ~ Max_life, correlation = corBrownian(phy = subset_tree),
                   data = test, method = "ML")
  

  output[count,] <- c(i, summary(pglsModel)$tTable[2,1],summary(pglsModel)$tTable[2,3],summary(pglsModel)$tTable[2,4])
}

pgls_results <- as_tibble(output) 
colnames(pgls_results) <- c("start","value","tvalue","pvalue")
pgls_results$scaf <- scaf
write_tsv(pgls_results, 
          paste0("../../data/gene_copies/aleutianus_genome/full_genome_scan_results/pgls_",scaf,".txt.gz"))
