library(tidyverse)
library("phyloseq")
library(ape)
library(phytools)
library(patchwork)
sample_order <- c("Sebastolobus_alascanus","Sebastes_iracundus","Sebastes_matsubarae","Sebastes_aleutianus",
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
                  "Sebastes_melanostomus","Sebastes_elongatus","Sebastes_saxicola","Sebastes_semicinctus","Sebastes_auriculatus","Sebastes_dallii",
                  "Sebastes_rastrelliger","Sebastes_nebulosus","Sebastes_atrovirens","Sebastes_aurora",
                  "Sebastes_caurinus","Sebastes_carnatus","Sebastes_maliger")
sample_order <- gsub("Sebastolobus","",gsub("Sebastes","",gsub("_","",sample_order)))


pdf("plots/sebonly_acti_per_gene_aligned_concat_distributions_20200806.pdf",height=10,width=4,useDingbats=FALSE)
dist_plot <- read_csv("mapping.txt") %>% 
  mutate(Map = case_when(Map == "X"~"0-0",
                         TRUE~ Map)) %>%
  separate(Map, c("start","end"),"-",convert=T) %>%
  ggplot(.,aes(x=start,xend=end,y=fct_rev(fct_relevel(Species, sample_order)),
               yend=fct_rev(fct_relevel(Species, sample_order)))) + 
  geom_segment(color="#5A9ECC",size=3) +
  theme_cowplot() +
  theme(axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank()) +
  labs(x = bquote('Distribution ('~x10^3~"km)"))
dist_plot
dev.off()

tree <- read.tree("../../data/trees/sebasto_sebaste_acti_per_gene_aligned_filt_ts_distfilt_genetrees/concat.treefile")
tree$tip.label <- gsub("Sebastes","S.",gsub("_"," ",tree$tip.label))
plot_tree(reroot(tree, 76, tree$edge.length[150]/2),ladderize=TRUE,label.tips="taxa_names", nodelabf=nodeplotboot())

pdf("plots/sebasto_sebaste_acti_per_gene_aligned_concat_tree_20200806.pdf",height=10,width=4,useDingbats=FALSE)
tree_plot <- plot_tree(reroot(tree, 76, 0.074),ladderize=TRUE,
                       label.tips="taxa_names", nodelabf=nodeplotboot()) +
  xlim(0,0.15)
tree_plot
dev.off()

#comparing different gene trees for all species or just sebastes
tree_seb <- read.tree("../../data/trees/sebasto_sebaste_acti_per_gene_aligned_filt_ts_distfilt_genetrees/concat.treefile")
tree_seb <- reroot(tree_seb, 76, 0.074)
tree_all <- read.tree("../../data/trees/all_species_acti_per_gene_aligned_filt_ts_distfilt_genetrees/concat.treefile")
both_species <- intersect(tree_seb$tip.label, tree_all$tip.label)
tree_all <- keep.tip(tree_all, both_species)
matched <- matrix(both_species,both_species)
double_species <- c(both_species,both_species)
dim(double_species) <- c(length(both_species), 2)
pdf("plots/phylogeny_comparison_20200806.pdf",height=20,width=20)
cophyloplot(tree_seb, tree_all,double_species,length.line=0,space=100,gap=20) +
  title("Sebastes vs All species")
dev.off()


