library(tidyverse)
library("phyloseq")
library(ape)
library(phytools)
library(patchwork)
sample_order <- c("Sebastiscus_albofasciatus","Sebastes_iracundus","Sebastes_matsubarae","Sebastes_baramenuke",
                  "Sebastes_aleutianus","Sebastes_glaucus","Sebastes_mentella","Sebastes_fasciatus","Sebastes_alutus",
                  "Sebastes_polyspinis","Sebastes_ciliatus","Sebastes_variabilis","Sebastes_itinus","Sebastes_minor",
                  "Sebastes_steindachneri","Sebastes_kiyomatsui","Sebastes_scythropus","Sebastes_thompsoni","Sebastes_inermis",
                  "Sebastes_joyneri","Sebastes_zonatus","Sebastes_trivittatus","Sebastes_taczanowskii","Sebastes_schlegelii",
                  "Sebastes_nivosus","Sebastes_hubbsi","Sebastes_oblongus","Sebastes_koreanus","Sebastes_pachycephalus",
                  "Sebastes_nudus","Sebastes_crameri","Sebastes_reedi","Sebastes_proriger","Sebastes_wilsoni",
                  "Sebastes_zacentrus","Sebastes_variegatus","Sebastes_melanops","Sebastes_flavidus","Sebastes_entomelas",
                  "Sebastes_mystinus","Sebastes_diaconus","Sebastes_ruberrimus","Sebastes_moseri","Sebastes_hopkinsi",
                  "Sebastes_melanostomus","Sebastes_elongatus","Sebastes_semicinctus","Sebastes_saxicola","Sebastes_auriculatus",
                  "Sebastes_rastrelliger","Sebastes_dallii","Sebastes_nebulosus","Sebastes_carnatus","Sebastes_atrovirens",
                  "Sebastes_maliger","Sebastes_aurora","Sebastes_caurinus","Sebastes_pinniger","Sebastes_miniatus",
                  "Sebastes_rosaceus","Sebastes_exsul","Sebastes_constellatus","Sebastes_oculatus","Sebastes_umbrosus",
                  "Sebastes_helvomaculatus","Sebastes_ensifer","Sebastes_rosenblatti","Sebastes_chlorostictus","Sebastes_levis",
                  "Sebastes_paucispinis","Sebastes_jordani","Sebastes_goodei","Sebastes_diploproa","Sebastes_babcocki",
                  "Sebastes_nigrocinctus","Sebastes_serriceps","Sebastes_rubrivinctus")
sample_order <- gsub("Sebastiscus","",gsub("Sebastes","",gsub("_","",sample_order)))


pdf("plots/sebonly_acti_per_gene_aligned_concat_distributions_20200421.pdf",height=10,width=4,useDingbats=FALSE)
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

tree <- read.tree("../../data/trees/sebonly_acti_per_gene_aligned_concat/sebonly_acti.renamed.treefile")
tree$tip.label <- gsub("Sebastes","S.",gsub("_"," ",tree$tip.label))
#tree <- midpoint.root(tree)
#tree <- root(tree, outgroup = "Sebastiscus albofasciatus",resolve.root = T)

pdf("plots/sebonly_acti_per_gene_aligned_concat_tree_20200421.pdf",height=10,width=4,useDingbats=FALSE)
tree_plot <- plot_tree(tree,ladderize=TRUE,label.tips="taxa_names", nodelabf=nodeplotboot()) +
  xlim(0,0.12)
tree_plot
dev.off()

tree <- phytools::reroot(tree,interactive=T) #Pick Sebastiscus albofasciatus
pdf("plots/sebonly_acti_per_gene_aligned_concat_tree_20200517.pdf",height=10,width=4,useDingbats=FALSE)
tree_plot <- plot_tree(tree,ladderize=TRUE,label.tips="taxa_names", nodelabf=nodeplotboot()) +
  xlim(0,0.15)
tree_plot
dev.off()

tree_plot + dist_plot 
