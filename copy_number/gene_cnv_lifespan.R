library(tidyverse)
library(cowplot)
library(qqman)
library(qvalue)
files <- list.files("../../data/gene_copies/aleutianus_genome/", "depth.gz")

gene_depths <- tibble(scaffold=character(),start=numeric(),end=numeric(),
                      gene=character(),total_depth=numeric(),
                      species=character(),sample=character())
for (i in files){
  species <- str_split(i,"\\.")[[1]][2]
  sample <- str_split(i,"\\.")[[1]][3]
  data <- read_tsv(paste0("../../data/gene_copies/aleutianus_genome/",i),
                   col_names = c("scaffold","start","end","gene","total_depth"))
  data$sample <- sample
  data$species <- species
  gene_depths <- rbind(gene_depths,data)
}
#Find top scaffolds
gene_depths %>%
  filter(sample == "A_latens_148487") %>%
  group_by(scaffold) %>%
  summarise(count=n()) %>%
  arrange(desc(count)) %>%
  head(24) %>% pull(scaffold) -> top_scaffolds

gene_depths %>%
  mutate(gene_length = end - start,
         average_depth = total_depth/gene_length) %>%
  group_by(sample) %>%
  summarize(median_depth = median(average_depth)) -> sample_depth
gene_depths %>%
  dplyr::select(species) %>%
  unique() %>%
  mutate(clade = case_when(grepl("Sebastes_", species) ~ "Sebastes",
                           TRUE ~ "Outgroup")) -> species_info
gene_depths %>%
  mutate(gene_length = end - start) %>%
  dplyr::select(gene,gene_length)  %>% 
  unique()-> gene_lengths

gene_depths %>%
  mutate(gene_length = end - start,
         average_depth = total_depth/gene_length) %>% 
  group_by(sample) %>%
  mutate(median_depth = median(average_depth)) %>%
  filter(species == "Sebastes_nigrocinctus") %>%
filter(average_depth > 100, average_depth < 1000) %>%
  ggplot(.,aes(x=average_depth)) + geom_density(fill="grey") +
  geom_vline(aes(xintercept=median_depth),linetype="dotted") + theme_cowplot()
  ungroup() %>%
  mutate(scaled_depth = (average_depth/median_depth)*2) %>%
  #filter(sample == "A_latens_148487") %>%
  #filter(scaffold %in% top_scaffolds) %>%
  mutate(binary_copy = case_when(scaled_depth < 0.5 ~ 0,
                                 scaled_depth < 1.5 ~ 1,
                                 scaled_depth <= 2.5 ~ 2,
                                 scaled_depth <= 3.5 ~ 3,
                                 scaled_depth <= 4.5 ~ 4,
                                 scaled_depth > 4.5 ~ 5)) %>%
  inner_join(species_info) %>%
  group_by(species,clade, sample,median_depth, binary_copy) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) ->summarized_copies
  
pdf("plots/CNV_sample_level_barchart_20200827.pdf",height=5,width=7)
summarized_copies %>%
  ggplot(aes(x=sample,y=freq, fill=as.factor(binary_copy))) + 
  geom_bar(position="stack", stat="identity") +
  theme_cowplot() +
  facet_wrap(~clade,scales="free_x") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_viridis_d(name="Copy number") +
  ylab("Frequency") + xlab("Sample")
dev.off()

pdf("plots/CNV_read_depth_20200827.pdf",height=5,width=7)
summarized_copies %>%
  filter(binary_copy == 2) %>%
  ggplot(.,aes(x=median_depth,y=freq)) +
  geom_point() +
  facet_wrap(~clade) +
  theme_cowplot() +
  ylab("Percent genes with 2 copies") +
  xlab("Median read depth")

summarized_copies %>%
  filter(binary_copy == 2) %>%
  ggplot(.,aes(x=median_depth,y=freq)) +
  geom_point() +
  facet_wrap(~clade) +
  theme_cowplot() +
  coord_cartesian(xlim=c(0,100)) +
  ylab("Percent genes with 2 copies") +
  xlab("Median read depth")
dev.off()

gene_depths %>%
    mutate(gene_length = end - start,
           average_depth = total_depth/gene_length) %>%
    group_by(sample) %>%
    mutate(median_depth = median(average_depth)) %>%
    ungroup() %>%
    mutate(scaled_depth = (average_depth/median_depth)*2) %>%
    mutate(binary_copy = case_when(scaled_depth < 0.5 ~ 0,
                                   scaled_depth < 1.5 ~ 1,
                                   scaled_depth <= 2.5 ~ 2,
                                   scaled_depth <= 3.5 ~ 3,
                                   scaled_depth <= 4.5 ~ 4,
                                   scaled_depth > 4.5 ~ 5)) %>%
  group_by(gene, binary_copy) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) ->summarized_genes

all_genes <- gene_depths$gene %>% unique()

summarized_genes %>%
  filter(binary_copy == 2) %>%
  arrange(freq) %>% pull(gene)-> summarized_gene_order

extra_genes <- all_genes[!all_genes %in% summarized_gene_order]

summarized_gene_order <- c(extra_genes, summarized_gene_order)

pdf("plots/CNV_gene_stats_20200827.pdf",height=5,width=7)
summarized_genes %>%
  ggplot(aes(fill=as.factor(binary_copy), y=freq, x=fct_relevel(gene, summarized_gene_order))) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis_d(name="Copy number") +
  xlab("Gene") +ylab("Frequency") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
  
summarized_genes %>%
  filter(binary_copy != 2) %>%
  group_by(gene) %>%
  summarize(total = sum(n),percent_not_diploid=total/97) %>%
  full_join(gene_lengths) %>%
  mutate(percent_not_diploid = case_when(is.na(percent_not_diploid) ~ 0,
                                         TRUE ~ percent_not_diploid)) %>%
  ggplot(.,aes(x=gene_length,y=percent_not_diploid)) +
  geom_point() +
  geom_smooth() +
  scale_x_continuous(trans='log10') +
  theme_cowplot()
  
summarized_genes %>%
  filter(binary_copy != 2) %>%
  group_by(gene) %>%
  summarize(total = sum(n),percent_not_diploid=total/97) %>%
  full_join(gene_lengths) %>%
  mutate(percent_not_diploid = case_when(is.na(percent_not_diploid) ~ 0,
                                         TRUE ~ percent_not_diploid)) %>%
  ggplot(.,aes(percent_not_diploid)) + geom_histogram() +
  theme_cowplot()
dev.off()


#####Using PGLS with nigrocinctus

max_life = read_tsv("../../data/meta/max_life_jul5_2020.txt") %>%
  filter(grepl("Sebastes_",species)) %>%
  filter(species != "Sebastes_aleutianus") %>%
  filter(!is.na(Max_life))

gene_depths %>%
  mutate(gene_length = end - start,
         average_depth = total_depth/gene_length) %>%
  group_by(sample) %>%
  mutate(median_depth = median(average_depth)) %>%
  ungroup() %>%
  mutate(scaled_depth = (average_depth/median_depth)*2) %>%
  mutate(binary_copy = case_when(scaled_depth <= 0.5 ~ 0,
                                 scaled_depth <= 1.5 ~ 1,
                                 scaled_depth <= 2.5 ~ 2,
                                 scaled_depth <= 3.5 ~ 3,
                                 scaled_depth <= 4.5 ~ 4,
                                 scaled_depth > 4.5 ~ 5)) %>%
  filter(species %in% max_life$species) %>%
  inner_join(max_life) -> processed_gene_depths
processed_gene_depths %>%
  group_by(gene, binary_copy) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  filter(binary_copy != 2) %>%
  group_by(gene) %>%
  summarize(total = sum(n),percent_not_diploid=total/97) %>%
  full_join(gene_lengths) %>%
  mutate(percent_not_diploid = case_when(is.na(percent_not_diploid) ~ 0,
                                         TRUE ~ percent_not_diploid)) %>%
  filter(percent_not_diploid > 0.05) %>%
  pull(gene)  -> Variable_genes

#Average value per species
gene_depths %>%
  mutate(gene_length = end - start,
         average_depth = total_depth/gene_length) %>%
  group_by(sample) %>%
  mutate(median_depth = median(average_depth)) %>%
  ungroup() %>%
  mutate(scaled_depth = (average_depth/median_depth)*2)  %>%
  group_by(gene,species) %>%
  summarize(ave_scaled_depth = mean(scaled_depth)) %>%
  inner_join(max_life) -> species_gene_depths

tree <- read.tree("../../data/trees/sebasto_sebaste_acti_per_gene_aligned_filt_ts_distfilt_genetrees/concat.treefile")

pgls_results_with_n <- tibble(gene=character(),value=numeric(),tvalue=numeric(),pvalue=numeric())

count<- 0
for (i in Variable_genes){
  count<- count+1
  print(count)
  species_gene_depths %>%
    filter(gene == i)  %>%
    filter(species %in% tree$tip.label) %>%
    as.data.frame -> test
  
  row.names(test) <- test$species
  pglsModel <- gls(ave_scaled_depth ~ Max_life, correlation = corBrownian(phy = keep.tip(tree, test$species)),
                   data = test, method = "ML")
  
  summary(pglsModel)$tTable[2,]
  tmp <- tibble(gene=i,value=summary(pglsModel)$tTable[2,1],
                tvalue=summary(pglsModel)$tTable[2,3],pvalue=summary(pglsModel)$tTable[2,4])
  pgls_results_with_n <- rbind(pgls_results_with_n, tmp)
}

pgls_results_with_n %>%
  ggplot(.,aes(pvalue)) + geom_histogram() +
  theme_cowplot()

#write_tsv(pgls_results_with_n, "../../data/gene_copies/sebastes_cnv_pgls.txt")
#write_tsv(species_gene_depths, "../../data/gene_copies/sebastes_cnv_pgls_inputdata.txt.gz")
species_gene_depths <- read_tsv("../../data/gene_copies/sebastes_cnv_pgls_inputdata.txt.gz")
pgls_results_with_n <- read_tsv("../../data/gene_copies/sebastes_cnv_pgls.txt")

pdf("plots/cnv_genes_qqplot.pdf",height=4,width=4,useDingbats = F)
qq(pgls_results_with_n$pvalue)
dev.off()
pgls_results_with_n$q <- qvalue(pgls_results_with_n$pvalue, fdr.level=0.2)$qvalues
#BH correction
pgls_results_with_n %>%
  arrange((pvalue)) %>%
  head(15) %>% pull(gene)-> top_genes

species_gene_depths %>%
  filter(gene %in% top_genes) %>%
  filter(species %in% tree$tip.label) %>% 
  ggplot(.,aes(x=Max_life,y=ave_scaled_depth)) + geom_point() +
  geom_smooth(method="lm") +
  theme_cowplot() +
  facet_wrap(~gene,scales="free_y")

pdf("plots/cnv_gene_butyrophilins.v1.pdf",height=3,width=6,useDingbats = F)
species_gene_depths %>%
  filter(gene == "FUN_060099" | gene == "FUN_060100") %>%
  filter(species %in% tree$tip.label) %>% 
  ggplot(.,aes(x=Max_life,y=ave_scaled_depth/2)) + geom_point(size=0.5) +
  geom_smooth(method="lm") +
  theme_cowplot() +
  facet_wrap(~gene,scales="free_y") +
  ylab("Normalized read depth") +
  theme(text=element_text(size=8),
        axis.text = element_text(size=8))
dev.off()
  

#####Using PGLS without nigrocinctus

max_life = read_tsv("../../data/meta/max_life_jul5_2020.txt") %>%
  filter(grepl("Sebastes_",species)) %>%
  filter(species != "Sebastes_aleutianus") %>%
  filter(species != "Sebastes_nigrocinctus") %>%
  filter(!is.na(Max_life))

gene_depths %>%
  mutate(gene_length = end - start,
         average_depth = total_depth/gene_length) %>%
  group_by(sample) %>%
  mutate(median_depth = median(average_depth)) %>%
  ungroup() %>%
  mutate(scaled_depth = (average_depth/median_depth)*2) %>%
  mutate(binary_copy = case_when(scaled_depth <= 0.5 ~ 0,
                                 scaled_depth <= 1.5 ~ 1,
                                 scaled_depth <= 2.5 ~ 2,
                                 scaled_depth <= 3.5 ~ 3,
                                 scaled_depth <= 4.5 ~ 4,
                                 scaled_depth > 4.5 ~ 5)) %>%
  filter(species %in% max_life$species) %>%
  inner_join(max_life) -> processed_gene_depths
processed_gene_depths %>%
  group_by(gene, binary_copy) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  filter(binary_copy != 2) %>%
  group_by(gene) %>%
  summarize(total = sum(n),percent_not_diploid=total/97) %>%
  full_join(gene_lengths) %>%
  mutate(percent_not_diploid = case_when(is.na(percent_not_diploid) ~ 0,
                                         TRUE ~ percent_not_diploid)) %>%
  filter(percent_not_diploid > 0.05) %>%
  pull(gene)  -> Variable_genes

#Average value per species
gene_depths %>%
  mutate(gene_length = end - start,
         average_depth = total_depth/gene_length) %>%
  group_by(sample) %>%
  mutate(median_depth = median(average_depth)) %>%
  ungroup() %>%
  mutate(scaled_depth = (average_depth/median_depth)*2)  %>%
  group_by(gene,species) %>%
  summarize(ave_scaled_depth = mean(scaled_depth)) %>%
  inner_join(max_life) -> species_gene_depths

tree <- read.tree("../../data/trees/sebasto_sebaste_acti_per_gene_aligned_filt_ts_distfilt_genetrees/concat.treefile")

pgls_results_without_n <- tibble(gene=character(),value=numeric(),tvalue=numeric(),pvalue=numeric())

count<- 0
for (i in Variable_genes){
  count<- count+1
  print(count)
  species_gene_depths %>%
    filter(gene == i)  %>%
    filter(species %in% tree$tip.label) %>%
    as.data.frame -> test
  
  row.names(test) <- test$species
  pglsModel <- gls(ave_scaled_depth ~ Max_life, correlation = corBrownian(phy = keep.tip(tree, test$species)),
                   data = test, method = "ML")
  
  summary(pglsModel)$tTable[2,]
  tmp <- tibble(gene=i,value=summary(pglsModel)$tTable[2,1],
                tvalue=summary(pglsModel)$tTable[2,3],pvalue=summary(pglsModel)$tTable[2,4])
  pgls_results_without_n <- rbind(pgls_results_without_n, tmp)
}

pgls_results_without_n %>%
  ggplot(.,aes(pvalue)) + geom_histogram() +
  theme_cowplot()

write_tsv(pgls_results_without_n, "../../data/gene_copies/sebastes_cnv_without_nigrocinctus_pgls.txt")
write_tsv(species_gene_depths, "../../data/gene_copies/sebastes_sebastes_cnv_without_nigrocinctus_pgls_inputdata.txt.gz")
qq(pgls_results_without_n$pvalue)
pgls_results_with_n$q <- qvalue(pgls_results_with_n$pvalue, pi0=1)$qvalues
#BH correction
pgls_results_with_n %>%
  arrange((q)) %>%
  head(16) %>% pull(gene)-> top_genes

species_gene_depths %>%
  filter(gene %in% top_genes) %>%
  filter(species %in% tree$tip.label) %>% 
  ggplot(.,aes(x=Max_life,y=ave_scaled_depth)) + geom_point() +
  geom_smooth(method="lm") +
  theme_cowplot() +
  facet_wrap(~gene,scales="free_y")


