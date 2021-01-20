#Plotting MSMC 

library(tidyverse)
library(ape)
library(nlme)
mutation_rate <- 2 *10^-9
#Mutation rate from "Moderate nucleotide diversity in the Atlantic herring is associated with a low mutation rate"
msmc <- tibble(time_index=character(),left_time_boundary=numeric(),right_time_boundary=numeric(),
               lambda_00=numeric(),sample=character(),species=character(),type=character(),assembly=character())
for (i in list.files("../../data/Ne/scaf_genomes/")){
  data <- read_tsv(paste0("../../data/Ne/scaf_genomes/",i))
  data$species <- str_split(i, "\\.")[[1]][1]
  data$sample <- str_split(i, "\\.")[[1]][2]
  data$type <- str_split(i, "\\.")[[1]][3]
  data$assembly <- "aleutianus_scaffolded"
  msmc <- rbind(msmc,data)
}
for (i in list.files("../../data/Ne/ref_genomes/")){
  data <- read_tsv(paste0("../../data/Ne/ref_genomes/",i))
  data$species <- str_split(i, "\\.")[[1]][1]
  data$sample <- str_split(i, "\\.")[[1]][2]
  data$type <- str_split(i, "\\.")[[1]][3]
  data$assembly <- "original"
  msmc <- rbind(msmc,data)
}
max_life = read_tsv("../../data/meta/max_life_jul5_2020.txt")


general_ranges <- read_tsv("../../data/meta/general_ranges.txt") %>%
  mutate(species = paste0(Genus,"_",Species))
write_tsv(msmc %>% mutate(type = case_when(type == "autosome" ~ "removed_2scaf",
                                           TRUE ~ "all_1MB_chr")), "../../data/Ne/allmsmc_2020_07_08.txt")
mutation_rate <- 2 *10^-9
msmc <- read_tsv("../../data/Ne/allmsmc_2020_07_08.txt")

pdf("plots/rockfish_refgenomes_Ne.v1.pdf",height=20,width=20)

msmc %>%
  group_by(sample,species) %>%
  filter(assembly=="aleutianus_scaffolded",type == "removed_2scaf") %>%
  filter(time_index > 1) %>%
  pivot_longer(c(-sample,-species,-lambda_00,-time_index,-type,-assembly), names_to = "boundary", values_to = "time_unscaled") %>%
  #filter(sample == "A_latens_148487") %>%
  mutate(generations=time_unscaled/mutation_rate,Ne=(1/lambda_00)/(2*mutation_rate),
         type.sample=paste0(type,".",sample)) %>%
  inner_join(.,general_ranges) %>%
  ggplot(aes(x=generations,y=Ne,color=species,group=type.sample,linetype=type)) + geom_line() +
  theme_cowplot() + scale_x_log10() +
  scale_y_log10() +
  ggtitle(paste0("Mu = ",mutation_rate)) +
  theme(legend.position = "none") +
  facet_wrap(~General_region) +
  coord_cartesian(ylim=c(5e4,1e7))

pdf("plots/rockfish_refgenomes_Ne.v2.pdf",height=30,width=30)
msmc %>%
  group_by(sample,species) %>%
  filter(time_index > 1) %>%
  filter(assembly=="aleutianus_scaffolded",type == "removed_2scaf") %>% 
  filter(species == "Sebastes_miniatus") %>% View()
  pivot_longer(c(-sample,-species,-lambda_00,-time_index,-type,-assembly), names_to = "boundary", values_to = "time_unscaled") %>%
  #filter(sample == "A_latens_148487") %>%
  mutate(generations=time_unscaled/mutation_rate,Ne=(1/lambda_00)/(2*mutation_rate),
         type.sample=paste0(type,".",sample)) %>%
  ggplot(aes(x=generations,y=Ne,color=species,group=type.sample,linetype=type)) + geom_line() +
  theme_cowplot() + scale_x_log10() +
  scale_y_log10() +
  ggtitle(paste0("Mu = ",mutation_rate)) +
  theme(legend.position = "none") +
  facet_wrap(~sample,scales="free_y")
dev.off()

#Find samples analyzed twice
repeated_samples <- msmc %>%
  filter(time_index == 0,type=="final") %>%
  group_by(sample) %>% count() %>% filter(n == 2) %>% pull(sample)

msmc %>%
  filter(sample %in% repeated_samples,type=="final") %>%
  group_by(sample,species,assembly) %>%
  pivot_longer(c(-sample,-species,-lambda_00,-time_index,-type,-assembly), names_to = "boundary", values_to = "time_unscaled") %>%
  mutate(generations=time_unscaled/mutation_rate,Ne=(1/lambda_00)/(2*mutation_rate),
         assembly.sample=paste0(assembly,".",sample)) %>%
  ggplot(aes(x=generations,y=Ne,color=species,group=assembly.sample,linetype=assembly)) + geom_line() +
  theme_cowplot() + scale_x_log10() +
  scale_y_log10() +
  ggtitle(paste0("Mu = ",mutation_rate)) +
  facet_wrap(~sample) +
  scale_x_log10() +
  scale_y_log10() 

dev.off()

msmc %>%
  filter(type =="removed_2scaf") %>%
  filter(species == "Sebastes_minor" | species == "Sebastes_variabilis" ) %>%
  filter(assembly == "aleutianus_scaffolded" ) %>%
  group_by(sample,species,assembly) %>%
  pivot_longer(c(-sample,-species,-lambda_00,-time_index,-type,-assembly), names_to = "boundary", values_to = "time_unscaled") %>%
  mutate(generations=time_unscaled/mutation_rate,Ne=(1/lambda_00)/(2*mutation_rate),
         assembly.sample=paste0(assembly,".",sample))  %>%
  ggplot(aes(x=generations,y=Ne,color=species,group=assembly.sample,linetype=assembly)) + geom_line() +
  theme_cowplot() + scale_x_log10() +
  scale_y_log10() +
  ggtitle(paste0("Mu = ",mutation_rate)) 

tree <- read.tree("../../data/trees/sebasto_sebaste_acti_per_gene_aligned_filt_ts_distfilt_genetrees/concat.treefile")


pgls_results_ne <- tibble(timing_i=numeric(),value=numeric(),tvalue=numeric(),pvalue=numeric())

for (i in 2:40){
  msmc %>%
    filter(type =="removed_2scaf") %>%
    filter(sample != "S-rosaceus_SEB-73", sample != "S-mystinus_SEB-17") %>%
    filter(time_index >=1, time_index < i) %>%
    filter(assembly == "aleutianus_scaffolded" ) %>%
    mutate(Ne=(1/lambda_00)/(2*mutation_rate)) %>%
    group_by(species,assembly) %>%
    summarize(mean_Ne  = mean(Ne)) %>% 
    inner_join(max_life) %>%
    filter(!is.na(Max_life)) %>%
    filter(grepl("Sebastes_",species)) %>%
    filter(species %in% tree$tip.label) %>%
    as.data.frame -> lifespan_ne
  
  lifespan_ne %>%  
    ggplot(.,aes(x=Max_life, y=mean_Ne)) +    geom_smooth(method="lm") +geom_point() +
    theme_cowplot() 
  
  
  
  row.names(lifespan_ne) <- lifespan_ne$species
  pglsModel <- gls(mean_Ne ~ Max_life, correlation = corBrownian(phy = keep.tip(tree, lifespan_ne$species)),
                   data = lifespan_ne, method = "ML")
  
  summary(pglsModel)$tTable[2,]
  tmp <- tibble(timing_i=i,value=summary(pglsModel)$tTable[2,1],
                tvalue=summary(pglsModel)$tTable[2,3],pvalue=summary(pglsModel)$tTable[2,4])
  pgls_results_ne <- rbind(pgls_results_ne, tmp)
}

pdf("plots/Lifespan_Ne_PGLS_range_20200825.pdf",height=3,width=5)
pgls_results_ne %>%
  ggplot(.,aes(x=timing_i,y=-log10(pvalue))) + geom_line() + theme_cowplot() +
  geom_hline(yintercept = -log10(0.05),linetype="dotted") +
  ylab("-log10(p-value)") +
  xlab("Timepoint limit")
dev.off()
  
i <- 20
msmc %>%
  filter(type =="removed_2scaf") %>%
  filter(sample != "S-rosaceus_SEB-73", sample != "S-mystinus_SEB-17") %>%
  filter(time_index >=1, time_index < i) %>%
  filter(assembly == "aleutianus_scaffolded" ) %>%
  mutate(Ne=(1/lambda_00)/(2*mutation_rate)) %>%
  group_by(species,assembly) %>%
  summarize(mean_Ne  = mean(Ne)) %>% 
  inner_join(max_life) %>%
  filter(!is.na(Max_life)) %>%
  filter(grepl("Sebastes_",species)) %>%
  filter(species %in% tree$tip.label) %>%
  as.data.frame -> lifespan_ne

write_tsv(lifespan_ne,"../../data/Ne/mean_Ne_i20_20201023.txt")
pvalue <- pgls_results_ne %>% filter(timing_i==20) %>% pull(pvalue)
pdf("plots/Lifespan_Ne_PGLS_20201023.pdf",height=3,width=5)
lifespan_ne %>%  
  ggplot(.,aes(x=Max_life, y=mean_Ne)) +    
  geom_smooth(method="lm",se=F,color="grey") +geom_point() +
  theme_cowplot() +
  ggtitle(paste0("PGLS p = ",formatC(pvalue, format = "e", digits = 2) )) +
  ylab("Mean Ne") + xlab("Max lifespan")

lifespan_ne %>%  
  mutate(quantile_rank = ntile(Max_life,4)) %>%
  ggplot(.,aes(x=as.factor(quantile_rank), y=mean_Ne)) +  
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width=0.1) +
  theme_cowplot() +
  ylab("Mean Ne") +
  xlab("Lifespan Quartile")
dev.off()

