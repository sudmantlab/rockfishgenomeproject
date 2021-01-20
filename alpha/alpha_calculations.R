library(tidyverse)
library(cowplot)
library(ape)
library(nlme)
library(patchwork)
#Load in all the het values

tree <- read.tree("../../data/trees/sebasto_sebaste_acti_per_gene_aligned_filt_ts_distfilt_genetrees/concat.treefile")

lifespan <- read_tsv("../../data/meta/max_life_jul5_2020.txt")
files <- list.files("../../data/alpha/het_annotation/")

het_data <- tibble(species=character(),sample=character(),gene=character(),
                   variants_effect_synonymous_variant=numeric(),
                   variants_effect_missense_variant=numeric())
for(file in files){
  tmp <- read_tsv(paste0("../../data/alpha/het_annotation/",file),skip=1) %>%
    rename(species=1,sample=2,gene=3) 
  if("variants_effect_synonymous_variant" %in% colnames(tmp)){
    if("variants_effect_missense_variant" %in% colnames(tmp)){
      tmp <- tmp %>% dplyr::select(species,sample,gene,variants_effect_synonymous_variant, variants_effect_missense_variant)
      het_data <- rbind(het_data,tmp)
    }

  }
  
}

#Load in all dn/ds tree data
gene_ids <- read_tsv("../../data/alpha/dnds_tree_values/seq_names.txt",col_names=c("orthogroup","gene")) %>%
  mutate(gene = gsub("-T1","",gene))

dnds_values <- read_tsv("../../data/alpha/dnds_tree_values/dnds_values.txt",
                        col_names=c("orthogroup","type","species","score")) %>%
  filter(score < 0.1) 

orthogroups_to_keep <- read_tsv("../../data/alpha/dnds_tree_values/dnds_usable_groups.txt",
                                col_names = "orthogroup")

ds_lengths <- read_tsv("../../data/alpha/dnds_tree_values/ds_lengths.txt",
                       col_names=c("orthogroup","length")) %>%
  mutate(type = "dS")
dn_lengths <- read_tsv("../../data/alpha/dnds_tree_values/dn_lengths.txt",
                       col_names=c("orthogroup","length")) %>%
  mutate(type = "dN")
dnds_values %>%
  inner_join(rbind(dn_lengths,ds_lengths)) %>%
  inner_join(orthogroups_to_keep) %>%
  mutate(count = length * score)%>%
  select(orthogroup,type,species,count)-> dnds_values
dnds_values <- gene_ids %>%
  inner_join(dnds_values) %>%
  pivot_wider(names_from=type,values_from=count) %>%
  filter(!is.na(dS)) %>%
  filter(!is.na(dN)) 

#Join them all together

dnds_values %>%
  inner_join(.,het_data) %>%
  mutate(alpha=1 - ((dS*variants_effect_missense_variant)/(dN*variants_effect_synonymous_variant)))-> all_data

pdf("plots/alpha_gene_counts.pdf",height=4,width=4)
all_data %>%
  group_by(orthogroup) %>%
  select(orthogroup,species) %>%
  summarize(n= n()) %>%
  ggplot(.,aes(n)) + geom_histogram() +
  theme_cowplot() + 
  ylab("Genes") +
  xlab("Species included") +
  geom_vline(xintercept=55.5,linetype="dashed") +
  geom_vline(xintercept=44,linetype="dashed")
dev.off()
all_data %>%
  group_by(orthogroup) %>%
  select(orthogroup,species) %>%
  summarize(n= n()) ->
  species_per_gene

filter(species_per_gene, n >= 44) -> called_genes

all_data %>%
  filter(orthogroup %in% called_genes$orthogroup) -> all_data
alpha_plot <- all_data %>% 
  group_by(species,sample) %>%
  summarize(dS=sum(dS),dN=sum(dN),syn_het=sum(variants_effect_synonymous_variant),
            non_het=sum(variants_effect_missense_variant),count=n()) %>%
  mutate(alpha=1 - ((dS*non_het)/(dN*syn_het))) %>% 
  inner_join(lifespan) %>%
  #filter(alpha > -3) %>%
  ggplot(.,aes(x=Max_life,y=alpha)) + geom_point() +
  theme_cowplot() +
  #coord_cartesian(ylim=c(-10,1)) +
  geom_hline(yintercept=1,linetype="dotted") 
  #geom_smooth(method="lm")

all_data %>% 
  group_by(species,sample) %>%
  summarize(dS=sum(dS),dN=sum(dN),syn_het=sum(variants_effect_synonymous_variant),
            non_het=sum(variants_effect_missense_variant),count=n()) %>%
  mutate(alpha=1 - ((dS*non_het)/(dN*syn_het))) %>% 
  group_by(species) %>%
  summarize(mean_alpha  = mean(alpha)) %>% 
  filter(species %in% tree$tip.label) %>%
  inner_join(lifespan) %>%
  filter(!is.na(Max_life)) %>%
  as.data.frame -> alpha_max_life
  
row.names(alpha_max_life) <- alpha_max_life$species
pglsModel <- gls(mean_alpha ~ Max_life, 
                 correlation = corBrownian(phy = keep.tip(tree, alpha_max_life$species),
                                           form = ~species),
                 data = alpha_max_life, method = "ML")
anova(pglsModel)

all_data %>% 
  group_by(species,sample) %>%
  summarize(dS=sum(dS),dN=sum(dN),syn_het=sum(variants_effect_synonymous_variant),
            non_het=sum(variants_effect_missense_variant),count=n()) %>%
  mutate(alpha=1 - ((dS*non_het)/(dN*syn_het))) %>% 
  inner_join(lifespan) %>%
  lm(alpha ~ Max_life, data=.) %>% summary()

all_data %>% 
  group_by(species,sample) %>%
  summarize(dS=sum(dS),dN=sum(dN),syn_het=sum(variants_effect_synonymous_variant),
            non_het=sum(variants_effect_missense_variant),count=n()) %>%
  mutate(dnds=dN/dS) %>% 
  inner_join(lifespan) %>% 
  ggplot(.,aes(x=Max_life,y=dnds)) + geom_point() +
  theme_cowplot() +
  #coord_cartesian(ylim=c(-10,1)) +
  geom_hline(yintercept=1,linetype="dotted") +
  ggtitle("Summed dN/dS")+
  geom_smooth(method="lm")

dnds_plot <- all_data %>% 
  group_by(species,sample) %>%
  summarize(dS=sum(dS),dN=sum(dN),syn_het=sum(variants_effect_synonymous_variant),
            non_het=sum(variants_effect_missense_variant),count=n()) %>%
  mutate(dnds=dN/dS) %>% 
  inner_join(lifespan) %>%
  ggplot(.,aes(x=Max_life,y=dnds)) + geom_point() +
  theme_cowplot() +
  ylab("Dn/Ds") +
  geom_smooth(method="lm",se=FALSE)
  



  
all_data %>% 
  group_by(species,sample) %>%
  summarize(dS=sum(dS),dN=sum(dN),syn_het=sum(variants_effect_synonymous_variant),
            non_het=sum(variants_effect_missense_variant),count=n()) %>%
  mutate(dnds=dN/dS) %>% 
  group_by(species) %>%
  summarize(mean_dnds  = mean(dnds)) %>% 
  filter(species %in% tree$tip.label) %>%
  inner_join(lifespan) %>%
  filter(!is.na(Max_life)) %>%
  as.data.frame -> dnds_max_life

row.names(dnds_max_life) <- dnds_max_life$species
pglsModel <- gls(mean_dnds ~ Max_life, 
                 correlation = corBrownian(phy = keep.tip(tree, dnds_max_life$species),
                                           form = ~species),
                 data = dnds_max_life, method = "ML")
summary(pglsModel)

pdf("plots/alpha_plots.pdf",height=4,width=6)
alpha_plot + dnds_plot + plot_annotation(tag_levels = 'A')
dev.off()

