library(ggtree)
library(treeio)
library(tidyverse)
library(ggtreeExtra)
tree <- treeio::read.beast("Desktop/allspecies_FigTree.20210311.tre")
tree@phylo$tip.label <- gsub("_"," ",tree@phylo$tip.label)
tree@phylo$tip.label <- gsub("Sebastes","S.",tree@phylo$tip.label)
tree@phylo$tip.label <- gsub("dalli","dallii",tree@phylo$tip.label)
tree@phylo$tip.label <- gsub("Hozukius","H.",tree@phylo$tip.label)
tree@phylo$tip.label <- gsub("Adelosebastes","A.",tree@phylo$tip.label)
tree@phylo$tip.label <- gsub("Sebastiscus","S.",tree@phylo$tip.label)
tree@phylo$tip.label <- gsub("Sebastolobus","S.",tree@phylo$tip.label)

lifespan <- read_tsv("Desktop/max_life_jul5_2020.txt") %>%
  mutate(shortened_life = case_when(is.na(Max_life) ~ 0,
                                    TRUE ~ Max_life/1)) %>%
  rename(ID=species) %>%
  mutate(ID = gsub("_", " ",ID)) %>%
  mutate(ID = gsub("Sebastes","S.",ID)) %>%
  mutate(ID = gsub("Hozukius","H.",ID)) %>%
  mutate(ID = gsub("Adelosebastes","A.",ID)) %>%
  mutate(ID = gsub("Sebastiscus","S.",ID)) %>%
  mutate(ID = gsub("Sebastolobus","S.",ID)) %>%
  mutate(ID = gsub("dalli","dallii",ID))


pdf("Documents/allspecies_FigTree.20210311.pdf",
    useDingbats = F,
    height=10,width=10)
ggtree(tree, layout="fan",open.angle=180)+ 
  geom_range(range='height_0.95_HPD', color='#6AADD5', size=2, alpha=.5) +
  geom_tiplab(size=1.5,hjust=-.1,
              align=T,
              linetype = "dotted") +
  geom_fruit(data=lifespan,
             geom=geom_bar,
             mapping=aes(y=ID, x=shortened_life,fill=Max_life),
             pwidth=1.3,
             stat="identity",
             orientation="y",
             position="auto",
             offset = 0.44) +
  scale_fill_viridis_c(option="plasma") 
dev.off()
