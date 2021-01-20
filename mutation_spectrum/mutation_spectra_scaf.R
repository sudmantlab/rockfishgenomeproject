library(tidyverse)
library(insect)
library(NMF)
library(broom)
library(nlme)
library(ape)
library(qvalue)
library(cowplot)
library(PNWColors)
library(caper)
colors <- pnw_palette("Bay",5)
file_list <- list.files(path="../../data/hetsites/scaf_genomes/",recursive=F)
max_life = read_tsv("../../data/meta/max_life_jul5_2020.txt")

species_order <- c("Sebastolobus_alascanus","Sebastolobus_altivelis",
                  "Sebastiscus_albofasciatus", "Sebastiscus_tertius","Adelosebastes_latens",
                  "Hozukius_emblemarius","Hozukius_guyotensis", "Sebastes_iracundus","Sebastes_matsubarae","Sebastes_aleutianus",
                  "Sebastes_baramenuke","Sebastes_alutus","Sebastes_polyspinis","Sebastes_ciliatus",
                  "Sebastes_variabilis","Sebastes_fasciatus","Sebastes_glaucus","Sebastes_mentella",
                  "Sebastes_itinus","Sebastes_minor","Sebastes_steindachneri",
                  "Sebastes_kiyomatsui","Sebastes_scythropus","Sebastes_thompsoni","Sebastes_inermis",
                  "Sebastes_joyneri","Sebastes_schlegelii","Sebastes_taczanowskii","Sebastes_trivittatus",
                  "Sebastes_zonatus","Sebastes_nivosus","Sebastes_hubbsi","Sebastes_oblongus","Sebastes_koreanus",
                  "Sebastes_nudus","Sebastes_pachycephalus","Sebastes_crameri","Sebastes_reedi",
                  "Sebastes_proriger","Sebastes_variegatus","Sebastes_wilsoni",
                  "Sebastes_zacentrus","Sebastes_flavidus",
                  "Sebastes_melanops","Sebastes_entomelas","Sebastes_diaconus",
                  "Sebastes_mystinus","Sebastes_miniatus","Sebastes_pinniger",
                  "Sebastes_ruberrimus","Sebastes_levis","Sebastes_jordani","Sebastes_goodei","Sebastes_paucispinis",
                  "Sebastes_babcocki","Sebastes_diploproa","Sebastes_nigrocinctus","Sebastes_rubrivinctus","Sebastes_serriceps",
                  "Sebastes_hopkinsi","Sebastes_moseri","Sebastes_rosaceus","Sebastes_constellatus","Sebastes_helvomaculatus",
                  "Sebastes_umbrosus","Sebastes_exsul","Sebastes_oculatus","Sebastes_ensifer","Sebastes_chlorostictus","Sebastes_rosenblatti",
                  "Sebastes_melanostomus","Sebastes_elongatus","Sebastes_saxicola","Sebastes_semicinctus","Sebastes_auriculatus","Sebastes_dalli",
                  "Sebastes_rastrelliger","Sebastes_nebulosus","Sebastes_atrovirens","Sebastes_aurora",
                  "Sebastes_caurinus","Sebastes_carnatus","Sebastes_maliger")

mutation_spectra <- tibble(mutation=character(),collapsed_n=numeric(),
                           species=character(),sample=character())
for (file in file_list){
  print(file)
  read_tsv(paste0("../../data/hetsites/scaf_genomes/", file)) %>%
    filter(!grepl("N",site_type )) %>%
      separate(site_type,c("ancestral","derived"),":") %>% 
      filter(!is.na(derived ))  %>%
      mutate(rev_ancestral = rc(ancestral),
             rev_derived = rc(derived),
             SNP_start = str_sub(ancestral,2,2)) %>%
      mutate(forward_form = paste0(ancestral,":",derived),
             rev_form = paste0(rev_ancestral,":",rev_derived)) %>%
      mutate(mutation = case_when((SNP_start == "A" | SNP_start == "C") ~ forward_form,
                                  TRUE ~ rev_form)) %>% 
      group_by(mutation) %>%
      summarize(collapsed_n = sum(n)) %>%
    mutate(filename=file) %>% separate(filename,c("species","sample"),"\\.") %>%
      select(mutation,collapsed_n,sample,species) -> tmp
    mutation_spectra <- rbind(mutation_spectra,tmp)
}


mutation_spectra$dataset <- "scaf"
#write_tsv(mutation_spectra,"../../data/hetsites/scaf_genomes/All_filtered.hetmutations.txt")

mutation_spectra <- read_tsv("../../data/hetsites/scaf_genomes/All_filtered.hetmutations.txt")
scaf_mutation_spectra <- mutation_spectra
pdf("plots/prelim_mutation_spectra.pdf",height=10,width=10)
mutation_spectra %>%
  group_by(sample,species) %>%
  mutate(total_count = sum(collapsed_n),
         proportion = collapsed_n/total_count)  %>%
  mutate(prime5 = str_sub(mutation,1,1),
         prime3 = str_sub(mutation,7,7),
         SNP_start = str_sub(mutation,2,2),
         SNP_end = str_sub(mutation,6,6),
         mut = paste0(SNP_start,"->",SNP_end))%>%
  ggplot(.,aes(y=paste0(mut,", ",prime5),x=prime3)) + geom_tile(aes(fill=proportion)) +
  scale_fill_viridis_c() +
  facet_wrap(~sample) +
  ylab("5-prime") + xlab("3-prime") + ggtitle("Prelim mutation spectra")

pdf("plots/mutation_spectra_for_phylogeny.pdf",height=20,width=20)
mutation_spectra %>%
  filter(species %in% species_order) %>%
  group_by(sample,species) %>%
  mutate(total_count = sum(collapsed_n),
         proportion = collapsed_n/total_count)  %>%
  mutate(prime5 = str_sub(mutation,1,1),
         prime3 = str_sub(mutation,7,7),
         SNP_start = str_sub(mutation,2,2),
         SNP_end = str_sub(mutation,6,6),
         mut = paste0(SNP_start,"->",SNP_end))%>%
  ungroup() %>%
  mutate(species = fct_relevel(species, species_order)) %>% 
  ggplot(.,aes(y=fct_reorder(sample, -as.numeric(species)),x=mutation)) + geom_tile(aes(fill=proportion)) +
  scale_fill_viridis_c() +
  facet_wrap(~mut, nrow=1,scales="free_x") +
  ylab("5-prime") + xlab("3-prime") + ggtitle("Mutation spectra") +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1))


mutation_spectra %>%
  filter(species %in% species_order) %>%
  group_by(sample,species) %>%
  mutate(total_count = sum(collapsed_n),
         proportion = collapsed_n/total_count)  %>%
  group_by(species,sample, mutation) %>%
  summarize(mean_prop = mean(proportion)) %>%
  group_by(mutation) %>%
  mutate(average_prop = mean(mean_prop), std_prop = stats::sd(mean_prop),z_prop = (mean_prop - average_prop)/std_prop)  %>%
  mutate(prime5 = str_sub(mutation,1,1),
         prime3 = str_sub(mutation,7,7),
         SNP_start = str_sub(mutation,2,2),
         SNP_end = str_sub(mutation,6,6),
         mut = paste0(SNP_start,"->",SNP_end),
         surround = paste0(prime5,"_",prime3))%>%
  inner_join(max_life) %>%
  ungroup() %>%
  mutate(species = fct_relevel(species, species_order)) %>% 
  ggplot(.,aes(y=fct_reorder(sample, -as.numeric(species)),x=mutation)) + geom_tile(aes(fill=z_prop)) +
  scale_fill_viridis_c() +
  facet_wrap(~mut, nrow=1,scales="free_x") +
  ylab("5-prime") + xlab("3-prime") + ggtitle("Mutation spectra normalized") +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1))

dev.off()


pdf("plots/mutation_spectra_all_types.pdf",height=7,width=7)
mutation_spectra %>%
  group_by(sample,species) %>%
  mutate(total_count = sum(collapsed_n),
         proportion = collapsed_n/total_count)  %>%
  mutate(prime5 = str_sub(mutation,1,1),
         prime3 = str_sub(mutation,7,7),
         SNP_start = str_sub(mutation,2,2),
         SNP_end = str_sub(mutation,6,6),
         mut = paste0(SNP_start,"->",SNP_end),
         surround = paste0(prime5,"_",prime3))%>%
  inner_join(max_life) %>%
  filter(!is.na(Max_life)) %>%
  filter(grepl("Sebastes",species)) %>%
  ggplot(.,aes(y=proportion,x=fct_reorder(surround,surround),
               group=fct_reorder(sample,Max_life),fill=Max_life)) + 
  geom_bar(position="dodge",stat="identity")+
  scale_fill_viridis_c() +
  facet_wrap(~mut,scales="free_x") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.2, hjust=1.2,
                                   size=7),
        legend.position="none") +
  xlab("Mutation") +
  ylab("Proportion") 
dev.off()
mutation_spectra %>%
  mutate(prime5 = str_sub(mutation,1,1),
         prime3 = str_sub(mutation,7,7),
         SNP_start = str_sub(mutation,2,2),
         SNP_end = str_sub(mutation,6,6),
         mut = paste0(SNP_start,"->",SNP_end),
         surround = paste0(prime5,"_",prime3))%>%
  group_by(sample,species,mut) %>%
  summarize(collapsed_n = sum(collapsed_n)) %>%
  group_by(sample,species) %>%
  mutate(total_count = sum(collapsed_n),
         proportion = collapsed_n/total_count)  %>%
  inner_join(max_life) %>%
  filter(!is.na(Max_life)) %>%
  ggplot(.,aes(y=proportion,x=mut,group=fct_reorder(sample,Max_life),fill=fct_reorder(species,Max_life))) + 
  geom_bar(position="dodge",stat="identity")+
  scale_fill_viridis_d() +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1),
        legend.position="none") +
  xlab("Mutation")




mutation_spectra %>%
  group_by(sample,species) %>%
  mutate(total_count = sum(collapsed_n),
         proportion = collapsed_n/total_count)  %>%
  group_by(species, mutation) %>%
  summarize(mean_prop = mean(proportion)) %>%
  group_by(mutation) %>%
  mutate(average_prop = mean(mean_prop), std_prop = stats::sd(mean_prop),z_prop = (mean_prop - average_prop)/std_prop)  %>%
  mutate(prime5 = str_sub(mutation,1,1),
         prime3 = str_sub(mutation,7,7),
         SNP_start = str_sub(mutation,2,2),
         SNP_end = str_sub(mutation,6,6),
         mut = paste0(SNP_start,"->",SNP_end),
         surround = paste0(prime5,"_",prime3))%>%
  inner_join(max_life) %>%
  #filter(!is.na(Max_life)) %>%
  ggplot(.,aes(y=fct_reorder(species,Max_life), x=fct_reorder(surround,surround),fill=z_prop)) +
  geom_tile() +
  facet_wrap(~mut,nrow=1) +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1),
        legend.position="none") 
  
mutation_spectra %>%
  mutate(prime5 = str_sub(mutation,1,1),
         prime3 = str_sub(mutation,7,7),
         SNP_start = str_sub(mutation,2,2),
         SNP_end = str_sub(mutation,6,6),
         mut = paste0(SNP_start,"->",SNP_end),
         surround = paste0(prime5,"_",prime3)) %>%
  group_by(sample, species, mut) %>%
  summarize(collapsed_n = sum(collapsed_n)) %>%
  group_by(sample,species) %>%
  mutate(total_count = sum(collapsed_n),
         proportion = collapsed_n/total_count)  %>%
  group_by(species, mut) %>%
  summarize(mean_prop = mean(proportion)) %>%
  group_by(mut) %>%
  mutate(average_prop = mean(mean_prop), std_prop = stats::sd(mean_prop),z_prop = (mean_prop - average_prop)/std_prop)  %>%
  inner_join(max_life) %>%
  filter(!is.na(Max_life)) %>%
  filter(!grepl("Sebastolobus",species)) %>%
  ggplot(.,aes(y=fct_reorder(species,Max_life), x=mut,fill=mean_prop)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1),
        legend.position="none") 

mutation_spectra %>%
  mutate(prime5 = str_sub(mutation,1,1),
         prime3 = str_sub(mutation,7,7),
         SNP_start = str_sub(mutation,2,2),
         SNP_end = str_sub(mutation,6,6),
         mut = paste0(SNP_start,"->",SNP_end),
         surround = paste0(prime5,"_",prime3)) %>%
  group_by(sample, species, mut) %>%
  summarize(collapsed_n = sum(collapsed_n)) %>%
  group_by(sample,species) %>%
  mutate(total_count = sum(collapsed_n),
         proportion = collapsed_n/total_count)  %>%
  group_by(species, mut) %>%
  summarize(mean_prop = mean(proportion)) %>%
  group_by(mut) %>%
  mutate(average_prop = mean(mean_prop), std_prop = stats::sd(mean_prop),z_prop = (mean_prop - average_prop)/std_prop)  %>%
  inner_join(max_life) %>%
  filter(!is.na(Max_life)) %>%
  filter(!grepl("Sebastolobus",species)) %>%
  ggplot(.,aes(x=Max_life, y=mean_prop)) +
  geom_point() + 
  geom_smooth(method="lm") +
  facet_wrap(~mut,scales="free_y") +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1),
        legend.position="none") +
  theme_cowplot()

mutation_spectra %>%
  mutate(prime5 = str_sub(mutation,1,1),
         prime3 = str_sub(mutation,7,7),
         SNP_start = str_sub(mutation,2,2),
         SNP_end = str_sub(mutation,6,6),
         mut = paste0(SNP_start,"->",SNP_end),
         surround = paste0(prime5,"_",prime3),
         last_two = paste0(mut,".", prime3)) %>%
  group_by(sample, species,last_two) %>%
  summarize(collapsed_n = sum(collapsed_n)) %>%
  group_by(sample,species) %>%
  mutate(total_count = sum(collapsed_n),
         proportion = collapsed_n/total_count)  %>%
  group_by(species,last_two) %>%
  summarize(mean_prop = mean(proportion)) %>%
  group_by(last_two) %>%
  mutate(average_prop = mean(mean_prop), std_prop = stats::sd(mean_prop),z_prop = (mean_prop - average_prop)/std_prop)  %>%
  inner_join(max_life) %>%
  filter(!is.na(Max_life)) %>%
  filter(!grepl("Sebastolobus",species)) %>%
  filter(last_two == "C->T.G") %>%
  ggplot(.,aes(x=Max_life, y=mean_prop)) +
  geom_point() + 
  geom_smooth(method="lm") +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1),
        legend.position="none") +
  theme_cowplot() +
  ggtitle(paste0("Proportion CpG mutations, p=",formatC(cpg_pvalue, format = "e", digits = 2)))


mutation_spectra %>%
  mutate(prime5 = str_sub(mutation,1,1),
         prime3 = str_sub(mutation,7,7),
         SNP_start = str_sub(mutation,2,2),
         SNP_end = str_sub(mutation,6,6),
         mut = paste0(SNP_start,"->",SNP_end),
         surround = paste0(prime5,"_",prime3),
         last_two = paste0(mut,".", prime3)) %>%
  group_by(sample, species,last_two) %>%
  summarize(collapsed_n = sum(collapsed_n)) %>%
  group_by(sample,species) %>%
  mutate(total_count = sum(collapsed_n),
         proportion = collapsed_n/total_count)  %>%
  group_by(species,last_two) %>%
  summarize(mean_prop = mean(proportion)) %>%
  group_by(last_two) %>%
  mutate(average_prop = mean(mean_prop), std_prop = stats::sd(mean_prop),z_prop = (mean_prop - average_prop)/std_prop)  %>%
  inner_join(max_life) %>%
  filter(!is.na(Max_life)) %>%
  filter(!grepl("Sebastolobus",species)) %>%
  do(tidy(lm(mean_prop~Max_life, .))) %>%
  filter(term == "Max_life") %>%
  filter(last_two == "C->T.G") %>% pull(p.value) -> cpg_pvalue
  
mutation_spectra %>%
  mutate(prime5 = str_sub(mutation,1,1),
         prime3 = str_sub(mutation,7,7),
         SNP_start = str_sub(mutation,2,2),
         SNP_end = str_sub(mutation,6,6),
         mut = paste0(SNP_start,"->",SNP_end),
         surround = paste0(prime5,"_",prime3)) %>%
  #group_by(sample, species, mut) %>%
  #summarize(collapsed_n = sum(collapsed_n)) %>%
  group_by(sample,species) %>%
  mutate(total_count = sum(collapsed_n),
         proportion = collapsed_n/total_count)  %>%
  group_by(species, mutation) %>%
  summarize(mean_prop = mean(proportion)) %>%
  group_by(mutation) %>%
  mutate(average_prop = mean(mean_prop), std_prop = stats::sd(mean_prop),z_prop = (mean_prop - average_prop)/std_prop)  %>%
  inner_join(max_life) %>%
  filter(!is.na(Max_life)) %>%
  filter(!grepl("Sebastolobus",species)) %>%
  do(tidy(lm(mean_prop~Max_life, .))) %>%
  filter(term == "Max_life") %>%   
  mutate(prime5 = str_sub(mutation,1,1),
                                          prime3 = str_sub(mutation,7,7),
                                          SNP_start = str_sub(mutation,2,2),
                                          SNP_end = str_sub(mutation,6,6),
                                          mut = paste0(SNP_start,"->",SNP_end),
                                          surround = paste0(prime5,"_",prime3)) %>% 
  dplyr::select(surround, mut, p.value) %>%
  arrange(mut) %>%
  mutate(permute="None")-> empirical_p


for (i in 1:1000){
  print(i)
  mutation_spectra %>%
    mutate(prime5 = str_sub(mutation,1,1),
           prime3 = str_sub(mutation,7,7),
           SNP_start = str_sub(mutation,2,2),
           SNP_end = str_sub(mutation,6,6),
           mut = paste0(SNP_start,"->",SNP_end),
           surround = paste0(prime5,"_",prime3)) %>%
    #group_by(sample, species, mut) %>%
    #summarize(collapsed_n = sum(collapsed_n)) %>%
    group_by(sample,species) %>%
    mutate(total_count = sum(collapsed_n),
           proportion = collapsed_n/total_count)  %>%
    group_by(species, mutation) %>%
    summarize(mean_prop = mean(proportion)) %>%
    group_by(mutation) %>%
    mutate(average_prop = mean(mean_prop), std_prop = stats::sd(mean_prop),z_prop = (mean_prop - average_prop)/std_prop)  %>%
    inner_join(max_life) %>%
    filter(!is.na(Max_life)) %>%
    filter(!grepl("Sebastolobus",species)) -> tmp
  tmp  %>%
    group_by(mutation) %>%
    mutate(Max_life=sample(Max_life)) %>%
    do(tidy(lm(mean_prop~Max_life, .))) %>%
    filter(term == "Max_life") %>%   
    mutate(prime5 = str_sub(mutation,1,1),
           prime3 = str_sub(mutation,7,7),
           SNP_start = str_sub(mutation,2,2),
           SNP_end = str_sub(mutation,6,6),
           mut = paste0(SNP_start,"->",SNP_end),
           surround = paste0(prime5,"_",prime3)) %>% 
    dplyr::select(surround, mut, p.value) %>%
    arrange(mut) %>% 
    mutate(permute = as.character(i)) -> permute_p
  empirical_p <- rbind(empirical_p, permute_p)
}
empirical_p %>%
  filter(permute != "None") %>%
  ggplot(.,aes(p.value)) +
  geom_density() +
  geom_vline(data=empirical_p %>% filter(permute == "None"),aes(xintercept=p.value)) +
  facet_wrap(~mutation) 

empirical_p %>%
  mutate(type = case_when(permute == "None" ~ "empirical",
                           TRUE ~ "permuted")) %>%
  ggplot(.,aes(p.value,fill=type)) +
  geom_density(alpha=0.2) +
  theme_cowplot() +
  scale_fill_manual(values=colors[c(1,5)])

#PCA
# ref_species <- c("Sebastes_entomelas","Sebastes_miniatus","Sebastes_pinniger","Sebastes_rosaceus",
#                  "Sebastes_schlegelii","Sebastes_umbrosus","Sebastolobus_alascanus")
# mutation_spectra %>%
#   filter(species %in% ref_species) -> mutation_spectra

mutation_spectra %>%
  group_by(sample,species) %>%
  mutate(total_count = sum(collapsed_n),
         proportion = collapsed_n/total_count) %>%
  ungroup() %>%
  select(-collapsed_n,-total_count,-species) %>%
  pivot_wider(names_from = sample, values_from = proportion) -> mutation_spectra_wide

mutation_spectra %>%
  select(sample,species) %>%
  unique() -> sample_info

mutation_spectra_pca <- prcomp(mutation_spectra_wide[c(3:ncol(mutation_spectra_wide))])

mutation_spectra_pca_results <-  as_tibble(mutation_spectra_pca$rotation) %>%
  mutate(sample = rownames(mutation_spectra_pca$rotation))

mutation_spectra_pca_loadings <-  as_tibble(mutation_spectra_pca$x) %>%
  mutate(mutation = mutation_spectra_wide$mutation)

mutation_spectra_pca_results %>%
  inner_join(sample_info) %>%
  ggplot(.,aes(x=PC1,y=PC2,color=species)) + 
  geom_point() +
  theme_cowplot()

mutation_spectra_pca_results %>%
  inner_join(sample_info) %>%
  inner_join(max_life) %>%
  lm(PC10 ~ Max_life, .) %>% summary()

mutation_spectra_pca_loadings %>%
  mutate(prime5 = str_sub(mutation,1,1),
         prime3 = str_sub(mutation,7,7),
         SNP_start = str_sub(mutation,2,2),
         SNP_end = str_sub(mutation,6,6),
         mut = paste0(SNP_start,"->",SNP_end)) %>%
  ggplot(.,aes(x=mutation,y=PC3)) + 
  geom_point() +
  facet_wrap(~mut,scales="free_x")


###Checking for different mutation associations with max age
#From "Rapid evolution of the human mutation spectrum"

all_mutations <- unique(mutation_spectra$mutation)


tree <- read.tree("../../data/trees/sebasto_sebaste_acti_per_gene_aligned_filt_ts_distfilt_genetrees/concat.treefile")
mutation_spectra %>%
  group_by(sample, species, mutation) %>%
  summarize(collapsed_n = sum(collapsed_n)) %>%
  group_by(sample,species) %>%
  mutate(total_count = sum(collapsed_n),
         proportion = collapsed_n/total_count)  %>%
  group_by(species, mutation) %>%
  summarize(species_prop = mean(proportion)) %>%
  group_by(mutation) %>%
  mutate(average_prop = mean(species_prop), std_prop = stats::sd(species_prop),z_prop = (species_prop - average_prop)/std_prop)  %>%
  inner_join(max_life) %>%
  filter(!is.na(Max_life)) %>%
  mutate(prime5 = str_sub(mutation,1,1),
         prime3 = str_sub(mutation,7,7),
         SNP_start = str_sub(mutation,2,2),
         SNP_end = str_sub(mutation,6,6),
         mut = paste0(SNP_start,"->",SNP_end),
         surround = paste0(prime5,"_",prime3)) %>%
  filter(!grepl("Sebastolobus",species)) -> mutation_spectra_species

mutation_spectra_species %>% pull (mutation) %>% unique() -> mutations
base_mutation_pvalue <- tibble(mutation=character(),value=numeric(),tvalue=numeric(),pvalue=numeric())
for (chosen_mutation in mutations){
  mutation_spectra_species %>%
    filter(mutation == chosen_mutation)  %>%
    filter(species %in% tree$tip.label) %>%
    as.data.frame -> test
  
  row.names(test) <- test$species
  pglsModel <- gls(z_prop ~ Max_life, correlation = corBrownian(phy = keep.tip(tree, test$species),
                                                                form = ~species),
                   data = test, method = "ML")
  
  summary(pglsModel)$tTable[2,]
  tmp <- tibble(mutation=chosen_mutation,value=summary(pglsModel)$tTable[2,1],
                tvalue=summary(pglsModel)$tTable[2,3],pvalue=summary(pglsModel)$tTable[2,4])
  base_mutation_pvalue <- rbind(base_mutation_pvalue, tmp)
}


base_mutation_pvalue <- base_mutation_pvalue %>% arrange(pvalue)


sequential_mutation_pvalue <- tibble(mutation = character(),pvalue = numeric(),value=numeric())

for (i in 2:(length(all_mutations)-1)){
  mutation_spectra %>%
    filter(!mutation %in% base_mutation_pvalue$mutation[1:(i-1)]) %>%
    group_by(sample, species, mutation) %>%
    summarize(collapsed_n = sum(collapsed_n)) %>%
    group_by(sample,species) %>%
    mutate(total_count = sum(collapsed_n),
           proportion = collapsed_n/total_count)  %>%
    group_by(species, mutation) %>%
    summarize(species_prop = mean(proportion)) %>%
    group_by(mutation) %>%
    mutate(average_prop = mean(species_prop), std_prop = stats::sd(species_prop),z_prop = (species_prop - average_prop)/std_prop)  %>%
    inner_join(max_life) %>%
    filter(!is.na(Max_life)) %>%
    mutate(prime5 = str_sub(mutation,1,1),
           prime3 = str_sub(mutation,7,7),
           SNP_start = str_sub(mutation,2,2),
           SNP_end = str_sub(mutation,6,6),
           mut = paste0(SNP_start,"->",SNP_end),
           surround = paste0(prime5,"_",prime3)) %>%
    filter(!grepl("Sebastolobus",species)) -> mutation_spectra_species_reduced
  mutation_spectra_species_reduced %>%
    filter(mutation == base_mutation_pvalue$mutation[i])  %>%
    filter(species %in% tree$tip.label) %>%
    as.data.frame -> test
  
  row.names(test) <- test$species
  pglsModel <- gls(z_prop ~ Max_life, correlation = corBrownian(phy = keep.tip(tree, test$species),
                                                                form = ~species),
                   data = test, method = "ML")
  
  summary(pglsModel)$tTable[2,]
  tmp <- tibble(mutation=base_mutation_pvalue$mutation[i],value=summary(pglsModel)$tTable[2,1],
                tvalue=summary(pglsModel)$tTable[2,3],pvalue=summary(pglsModel)$tTable[2,4])

  sequential_mutation_pvalue <- rbind(sequential_mutation_pvalue,tmp)
}
sequential_mutation_pvalue <- sequential_mutation_pvalue %>% rbind(base_mutation_pvalue[1,],.) %>%
  rbind(.,base_mutation_pvalue[96,])


sequential_mutation_pvalue$qvalue <- qvalue(sequential_mutation_pvalue$pvalue)$qvalues

sequential_mutation_pvalue %>%
  mutate(significant = case_when(qvalue < 0.05 ~ "sig",
                                 TRUE ~ "nonsig")) %>%
  mutate(prime5 = str_sub(mutation,1,1),
         prime3 = str_sub(mutation,7,7),
         SNP_start = str_sub(mutation,2,2),
         SNP_end = str_sub(mutation,6,6),
         mut = paste0(SNP_start,"->",SNP_end)) -> tmp_plot_data

pdf("plots/mutation_spectra_lifemax_sequential_pgls_v2.pdf",height=5,width=3,useDingbats=FALSE)
tmp_plot_data %>%
  ggplot(.,aes(y=paste0(mut,", ",prime5),x=prime3)) + geom_tile(aes(fill=value)) +
  scale_fill_distiller(palette = "Spectral",name="Enrichment") +
  ylab("5-prime") + xlab("3-prime") + 
  #ggtitle("Effect of max lifespan on mutation type") +
  geom_point(data = tmp_plot_data %>% filter(significant == "sig"),
             shape=8) +
  theme_cowplot()
dev.off()

pdf("plots/mutation_spectra_lifemax_sequential_pgls_v2L.pdf",height=2,width=7,useDingbats=FALSE)
tmp_plot_data %>%
  ggplot(.,aes(x=paste0(mut,", ",prime5),y=prime3)) + geom_tile(aes(fill=value)) +
  scale_fill_distiller(palette = "Spectral",name="Enrichment") +
  ylab("5-prime") + xlab("3-prime") + 
  #ggtitle("Effect of max lifespan on mutation type") +
  geom_point(data = tmp_plot_data %>% filter(significant == "sig"),
             shape=8) +
  theme_cowplot() +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  theme(axis.text = element_text(size=7))
  
dev.off()


####PGLS

tree <- read.tree("../../data/trees/sebasto_sebaste_acti_per_gene_aligned_filt_ts_distfilt_genetrees/concat.treefile")
mutation_spectra %>%
  group_by(sample, species, mutation) %>%
  summarize(collapsed_n = sum(collapsed_n)) %>%
  group_by(sample,species) %>%
  mutate(total_count = sum(collapsed_n),
         proportion = collapsed_n/total_count)  %>%
  group_by(species, mutation) %>%
  summarize(species_prop = mean(proportion)) %>%
  group_by(mutation) %>%
  mutate(average_prop = mean(species_prop), std_prop = stats::sd(species_prop),z_prop = (species_prop - average_prop)/std_prop)  %>%
  inner_join(max_life) %>%
  filter(!is.na(Max_life)) %>%
  mutate(prime5 = str_sub(mutation,1,1),
         prime3 = str_sub(mutation,7,7),
         SNP_start = str_sub(mutation,2,2),
         SNP_end = str_sub(mutation,6,6),
         mut = paste0(SNP_start,"->",SNP_end),
         surround = paste0(prime5,"_",prime3)) %>%
  filter(!grepl("Sebastolobus",species)) -> mutation_spectra_species

mutation_spectra_species %>% pull (mutation) %>% unique() -> mutations
pgls_results <- tibble(mutation=character(),value=numeric(),tvalue=numeric(),pvalue=numeric())
for (chosen_mutation in mutations){
  mutation_spectra_species %>%
    filter(mutation == chosen_mutation)  %>%
    filter(species %in% tree$tip.label) %>%
    as.data.frame -> test
  
  row.names(test) <- test$species
  pglsModel <- gls(z_prop ~ Max_life, correlation = corBrownian(phy = keep.tip(tree, test$species),
                                                                form = ~species),
                   data = test, method = "ML")
  
  summary(pglsModel)$tTable[2,]
  tmp <- tibble(mutation=chosen_mutation,value=summary(pglsModel)$tTable[2,1],
                tvalue=summary(pglsModel)$tTable[2,3],pvalue=summary(pglsModel)$tTable[2,4])
  pgls_results <- rbind(pgls_results, tmp)
}

pgls_results$qvalue <- qvalue(pgls_results$pvalue)$qvalues

pgls_results %>%
  mutate(significant = case_when(qvalue < 0.05 ~ "sig",
                                 TRUE ~ "nonsig")) %>%
  mutate(prime5 = str_sub(mutation,1,1),
         prime3 = str_sub(mutation,7,7),
         SNP_start = str_sub(mutation,2,2),
         SNP_end = str_sub(mutation,6,6),
         mut = paste0(SNP_start,"->",SNP_end)) -> tmp_plot_data
  
pdf("plots/mutation_spectra_lifemax_pgls_v2.pdf",height=5,width=3,useDingbats=FALSE)
tmp_plot_data %>%
  ggplot(.,aes(y=paste0(mut,", ",prime5),x=prime3)) + geom_tile(aes(fill=value)) +
  scale_fill_distiller(palette = "Spectral") +
  ylab("5-prime") + xlab("3-prime") + 
  #ggtitle("Effect of max lifespan on mutation type") +
  geom_point(data = tmp_plot_data %>% filter(significant == "sig"),
            shape=8) +
  theme_cowplot()
dev.off()

pdf("plots/mutation_spectra_lifemax_pgls_v3L.pdf",height=1.6,width=7,useDingbats=FALSE)
tmp_plot_data %>%
  ggplot(.,aes(x=paste0(mut,", ",prime5),y=prime3)) + geom_tile(aes(fill=value)) +
  scale_fill_distiller(palette = "Spectral") +
  ylab("5-prime") + xlab("3-prime") + 
  #ggtitle("Effect of max lifespan on mutation type") +
  geom_point(data = tmp_plot_data %>% filter(significant == "sig"),
             shape=8) +
  theme_cowplot() +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  theme(text = element_text(size = 7),
        axis.text = element_text(size=7)) 
dev.off()
  
#PGLS for CpG
mutation_spectra %>%
  mutate(prime5 = str_sub(mutation,1,1),
         prime3 = str_sub(mutation,7,7),
         SNP_start = str_sub(mutation,2,2),
         SNP_end = str_sub(mutation,6,6),
         mut = paste0(SNP_start,"->",SNP_end),
         surround = paste0(prime5,"_",prime3),
         last_two = paste0(mut,".", prime3)) %>%
  group_by(sample, species,last_two) %>%
  summarize(collapsed_n = sum(collapsed_n)) %>%
  group_by(sample,species) %>%
  mutate(total_count = sum(collapsed_n),
         proportion = collapsed_n/total_count)  %>%
  group_by(species,last_two) %>%
  summarize(mean_prop = mean(proportion)) %>%
  group_by(last_two) %>%
  mutate(average_prop = mean(mean_prop), std_prop = stats::sd(mean_prop),z_prop = (mean_prop - average_prop)/std_prop)  %>%
  inner_join(max_life) %>%
  filter(!is.na(Max_life)) %>%
  filter(species %in% tree$tip.label) %>%
  filter(last_two == "C->T.G") %>%
  as.data.frame -> test

row.names(test) <- test$species
pglsModel <- gls(mean_prop ~ Max_life, correlation = corBrownian(phy = keep.tip(tree, test$species)),
                 data = test, method = "ML")

summary(pglsModel)$tTable[2,4]
test %>% 
  filter(grepl("Sebastes_",species)) %>%
ggplot(.,aes(x=Max_life, y=mean_prop)) +
  geom_point() + 
  geom_abline(slope=coef(pglsModel)[2],intercept=coef(pglsModel)[1]) +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1),
        legend.position="none") +
  theme_cowplot() +
  ylab("CpG mutation proportion") +
  ggtitle(paste0("Proportion CpG mutations\nPGLS p=",formatC(summary(pglsModel)$tTable[2,4], format = "e", digits = 2)))

      
#LM for CpG

mutation_spectra %>%
  mutate(prime5 = str_sub(mutation,1,1),
         prime3 = str_sub(mutation,7,7),
         SNP_start = str_sub(mutation,2,2),
         SNP_end = str_sub(mutation,6,6),
         mut = paste0(SNP_start,"->",SNP_end),
         surround = paste0(prime5,"_",prime3),
         last_two = paste0(mut,".", prime3)) %>%
  group_by(sample, species,last_two) %>%
  summarize(collapsed_n = sum(collapsed_n)) %>%
  group_by(sample,species) %>%
  mutate(total_count = sum(collapsed_n),
         proportion = collapsed_n/total_count)  %>%
  group_by(species,last_two) %>%
  summarize(mean_prop = mean(proportion)) %>%
  group_by(last_two) %>%
  mutate(average_prop = mean(mean_prop), std_prop = stats::sd(mean_prop),z_prop = (mean_prop - average_prop)/std_prop)  %>%
  inner_join(max_life) %>%
  filter(!is.na(Max_life)) %>%
  filter(!grepl("Sebastolobus",species)) %>%
  do(tidy(lm(mean_prop~Max_life, .))) %>%
  filter(term == "Max_life") %>%
  filter(last_two == "C->T.G") %>% pull(p.value) -> cpg_pvalue

mutation_spectra %>%
  mutate(prime5 = str_sub(mutation,1,1),
         prime3 = str_sub(mutation,7,7),
         SNP_start = str_sub(mutation,2,2),
         SNP_end = str_sub(mutation,6,6),
         mut = paste0(SNP_start,"->",SNP_end),
         surround = paste0(prime5,"_",prime3),
         last_two = paste0(mut,".", prime3)) %>%
  group_by(sample, species,last_two) %>%
  summarize(collapsed_n = sum(collapsed_n)) %>%
  group_by(sample,species) %>%
  mutate(total_count = sum(collapsed_n),
         proportion = collapsed_n/total_count)  %>%
  group_by(species,last_two) %>%
  summarize(mean_prop = mean(proportion)) %>%
  group_by(last_two) %>%
  mutate(average_prop = mean(mean_prop), std_prop = stats::sd(mean_prop),z_prop = (mean_prop - average_prop)/std_prop)  %>%
  inner_join(max_life) %>%
  filter(!is.na(Max_life)) %>%
  filter(!grepl("Sebastolobus",species)) %>%
  filter(last_two == "C->T.G") -> plotting_data

pdf("plots/mutation_spectra_lifemax_pgls_cpg_v1.pdf",height=3,width=5)
plotting_data %>%
  ggplot(.,aes(x=Max_life, y=mean_prop)) +
  geom_point() + 
  geom_smooth(method="lm",se=F,color="black") +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1),
        legend.position="none") +
  theme_cowplot() +
  ylab("CpG mutation proportion") +
  xlab("Maximum lifespan") +
  # geom_segment(aes(x=min(plotting_data$Max_life),xend=max(plotting_data$Max_life),
  #                  y=(min(plotting_data$Max_life)*coef(pglsModel)[2] + coef(pglsModel)[1]),
  #                  yend=(max(plotting_data$Max_life)*coef(pglsModel)[2] + coef(pglsModel)[1])),
  #              linetype="dotted") +
  ggtitle(paste0("LM p=",formatC(cpg_pvalue, format = "e", digits = 2),"\nPGLS p=",formatC(summary(pglsModel)$tTable[2,4], format = "e", digits = 2)))
dev.off()



  
         