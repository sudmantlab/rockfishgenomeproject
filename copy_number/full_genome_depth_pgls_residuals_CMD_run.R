
library(tidyverse)
library(cowplot)
library(nlme)
library(ape)
args = commandArgs(trailingOnly=TRUE)

scaf <- args[1]
scaf <- "PGA_scaffold_74__15"
sample_depth <- read_tsv("../../data/gene_copies/average_sample_depth.txt")

residuals = read_tsv("../../data/meta/Lmat_50_depth_model_49species.txt") %>%
  mutate(species = gsub("S. ","Sebastes_",species))

tree <- read.tree("../../data/trees/sebasto_sebaste_acti_per_gene_aligned_filt_ts_distfilt_genetrees/concat.treefile")




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
  filter(start > 37000000, start < 39000000) %>%
  inner_join(.,sample_depth) %>%
  filter(species  != "Sebastes_aleutianus") %>%
  mutate(scaled_depth = total_depth/median_depth/100) %>%
  group_by(scaffold,start,end,species) %>%
  summarize(mean_scaled_depth = mean(scaled_depth)) %>%
  ungroup() %>%
  dplyr::select(scaffold, start, end,species, mean_scaled_depth) %>%
  inner_join(residuals) %>%
  filter(species %in% tree$tip.label) -> region_depth_selected

species_used <- unique(region_depth_scaled$species)
subset_tree <- keep.tip(tree, species_used)



output <- matrix(ncol=10, nrow=length(unique(region_depth_selected$start)))
count<- 0
for (i in unique(region_depth_selected$start)){
#for (i in head(unique(region_depth_scaled$start),100)){
  
  count<- count+1
  print(count)
  region_depth_selected %>%
    filter(start == i)  %>%
    as.data.frame -> test
  if (var(test$mean_scaled_depth) == 0){next;}
  
  row.names(test) <- test$species
  pglsModel_depth <- gls(mean_scaled_depth ~ max_depth, correlation = corBrownian(phy = subset_tree),
                   data = test, method = "ML")
  pglsModel_resid <- gls(mean_scaled_depth ~ m_depth_matL_residual, correlation = corBrownian(phy = subset_tree),
                         data = test, method = "ML")
  pglsModel_size <- gls(mean_scaled_depth ~ Lmat50, correlation = corBrownian(phy = subset_tree),
                         data = test, method = "ML")
  

  output[count,] <- c(i, summary(pglsModel_depth)$tTable[2,1],summary(pglsModel_depth)$tTable[2,3],summary(pglsModel_depth)$tTable[2,4],
                      summary(pglsModel_size)$tTable[2,1],summary(pglsModel_size)$tTable[2,3],summary(pglsModel_size)$tTable[2,4],
                      summary(pglsModel_resid)$tTable[2,1],summary(pglsModel_resid)$tTable[2,3],summary(pglsModel_resid)$tTable[2,4])
}

pgls_results <- as_tibble(output) 
colnames(pgls_results) <- c("start","value_depth","tvalue_depth","pvalue_depth",
                            "value_size","tvalue_size","pvalue_size",
                            "value_resid","tvalue_resid","pvalue_resid")
pgls_results$scaf <- scaf
write_tsv(pgls_results, 
          paste0("../../data/gene_copies/aleutianus_genome/full_genome_scan_results/pgls_resid_",scaf,".txt.gz"))
