library(tidyverse)
library(zoo)
library(PNWColors)
library(cowplot)

colors <- pnw_palette("Bay",5)
files <- list.files("../../data/gene_copies/aleutianus_genome/full_genome_scan_results/")
chr_names <- read_tsv("../../data/meta/chromosome_numbers.txt") %>%
  rename(chr_n = X7) %>% 
  mutate(scaf = gsub("Sebastes_aleutianus.","",CHR)) %>%
  select(scaf,chr_n) 
all_data <- tibble(start=character(),value=numeric(),tvalue=numeric(),pvalue=numeric(),scaf=character())
for (i in files){
  data <- read_tsv(paste0("../../data/gene_copies/aleutianus_genome/full_genome_scan_results/",i))
  all_data <- rbind(all_data,data)
}

all_data <- all_data %>%
  inner_join(chr_names) %>%
  mutate(log_pvalue = case_when(value > 0 ~ -log10(pvalue),
                                TRUE ~ log10(pvalue))) 

chr_lengths <- all_data %>%
  group_by(chr_n) %>%
  summarise(length=max(start,na.rm=T)) 

all_data_cum <- chr_lengths %>% 
  select(chr_n,length) %>%
  rename(chr_len = length) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(all_data, ., by=c("chr_n"="chr_n")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr_n, start) %>%
  mutate( start_cum=start+tot)

axisdf = all_data_cum %>% group_by(chr_n) %>% summarize(center=( max(start_cum,na.rm=T) + min(start_cum,na.rm=T) ) / 2 )



but_genes <- read_tsv("gene_coords_all-butro.gff3.txt",
                      col_names = c("scaf_full","V1","V2","begin","end","V3","V4","V5","V6")) %>%
  mutate(scaf = gsub("Sebastes_aleutianus.","",scaf_full)) %>%
  select(scaf,begin,end) %>%
  inner_join(chr_names) %>%
  filter(chr_n %in% unique(all_data$chr_n)) %>%
  mutate(mid = (begin+end)/2) %>%
  select(chr_n,mid)

but_genes <- chr_lengths %>% 
  select(chr_n,length) %>%
  rename(chr_len = length) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(but_genes, ., by=c("chr_n"="chr_n")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr_n, mid) %>%
  mutate( mid_cum=mid+tot)



p_threshold = 3

pdf("plots/read_depth_pgls_gwas.v1.pdf",height=3.5,width=15)

all_data_cum %>%
  mutate(window = floor(start_cum/10000)*10000) %>%
  group_by(chr_n,window) %>%
  summarize(mean_log_p = median(abs(log_pvalue)),
            high_count = sum(log_pvalue > p_threshold)/n()) %>%
  ggplot(.) + 
  geom_line( aes(x=window, y=high_count,color=as.factor(chr_n)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c(colors[1], colors[3]), 24 )) +
  theme_cowplot() +
  ylab(paste0("Proportion regions log(pvalue) > ",p_threshold)) +
  ggtitle("Positive") +
  theme(legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank()) +
  scale_x_continuous(label = axisdf$chr_n, 
                      breaks= axisdf$center,
                      guide = guide_axis(n.dodge=2) ) +
  coord_cartesian(ylim=c(0,0.5)) +
  geom_point(data=but_genes,aes(x=mid_cum,y=0),color="red",shape=124)


all_data_cum %>%
  mutate(window = floor(start_cum/10000)*10000) %>%
  group_by(chr_n,window) %>%
  summarize(mean_log_p = median(abs(log_pvalue)),
            high_count = sum(log_pvalue < -p_threshold)/n()) %>%
  ggplot(.) + 
  geom_line( aes(x=window, y=high_count,color=as.factor(chr_n)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c(colors[1], colors[3]), 24 )) +
  theme_cowplot() +
  ylab(paste0("Proportion regions log(pvalue) > ",p_threshold)) +
  ggtitle("Negative") +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank()) +
  scale_x_continuous( label = axisdf$chr_n, breaks= axisdf$center,guide = guide_axis(n.dodge=2) ) +
  coord_cartesian(ylim=c(0,0.5)) +
  geom_point(data=but_genes,aes(x=mid_cum,y=0),color="red",shape=124)


all_data_cum %>%
  mutate(window = floor(start_cum/10000)*10000) %>%
  group_by(chr_n,window) %>%
  summarize(mean_log_p = median(abs(log_pvalue)),
            high_count = sum(log_pvalue > p_threshold)/n()) -> positive_windows

all_data_cum %>%
  mutate(window = floor(start_cum/10000)*10000) %>%
  group_by(chr_n,window) %>%
  summarize(mean_log_p = median(abs(log_pvalue)),
            high_count = -sum(log_pvalue < -p_threshold)/n()) -> negative_windows

positive_windows %>%
ggplot(.) + 
  geom_line( aes(x=window, y=high_count,color=as.factor(chr_n)), alpha=0.8, size=1) +
  geom_line(data=negative_windows, aes(x=window, y=high_count,color=as.factor(chr_n)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c(colors[1], colors[3]), 24 )) +
  theme_cowplot() +
  ylab(paste0("Proportion regions log(pvalue) > ",p_threshold)) +
  theme( legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank()) +
  scale_x_continuous( label = axisdf$chr_n, breaks= axisdf$center,guide = guide_axis(n.dodge=2) ) +
  coord_cartesian(ylim=c(-0.1,0.5)) +
  geom_point(data=but_genes,aes(x=mid_cum,y=0),color="red",shape=124) 

all_data_cum %>%
  filter(2.5 <= log_pvalue) %>%
  ggplot(.) +
  
  # Show all points
  geom_point( aes(x=start_cum, y=log_pvalue,color=as.factor(chr_n)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c(colors[1], colors[3]), 24 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr_n, breaks= axisdf$center,guide = guide_axis(n.dodge=2)  ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw() +
  theme( legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank()) +
  coord_cartesian(ylim=c(2.5,8)) +
  ggtitle("Positive" ) +
  xlab("Chromosome") +
  ylab("Log pvalue") +
  geom_point(data=but_genes,aes(x=mid_cum,y=2.5),color="red",shape=124)

all_data_cum %>%
  filter(0 <= log_pvalue) %>%
  ggplot(.) +
  
  # Show all points
  geom_point( aes(x=start_cum, y=log_pvalue,color=as.factor(chr_n)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c(colors[1], colors[3]), 24 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr_n, breaks= axisdf$center,guide = guide_axis(n.dodge=2)  ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw() +
  theme( legend.position="none",
         panel.border = element_blank(),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         axis.title.x = element_blank()) +
  ggtitle("Positive" ) +
  xlab("Chromosome") +
  ylab("Log pvalue") +
  geom_point(data=but_genes,aes(x=mid_cum,y=0),color="red",shape=124) +
  coord_cartesian(xlim=c(37282801, 38282801))



all_data_cum %>%
  filter(-2.5 >= log_pvalue) %>%
  ggplot(.) +
  
  # Show all points
  geom_point( aes(x=start_cum, y=abs(log_pvalue),color=as.factor(chr_n)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c(colors[1], colors[3]), 24 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr_n, breaks= axisdf$center,guide = guide_axis(n.dodge=2)  ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw() +
  theme( legend.position="none",
         panel.border = element_blank(),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         axis.title.x = element_blank()) +
  coord_cartesian(ylim=c(2.5,8)) +
  ggtitle("Negative" ) +
  xlab("Chromosome") +
  ylab("Log pvalue") +
  geom_point(data=but_genes,aes(x=mid_cum,y=2.5),color="red",shape=124)




all_data %>%
  mutate(direction = case_when(tvalue < 0 ~ "Negative",
                               TRUE ~ "Positive"),
         pvalue_bin = floor(abs(log_pvalue)*10)/10) %>%
  filter(pvalue_bin >= 3) %>%
  ggplot(.,aes(x=pvalue_bin,fill=direction)) + 
  geom_bar(position="dodge2",stat="count") +
  theme_cowplot() +
  xlab("Log P-value bin") +
  scale_fill_manual(values = c(colors[2], colors[5]),name="Correlation direction")
  
all_data %>%
  mutate(direction = case_when(tvalue < 0 ~ "Negative",
                               TRUE ~ "Positive"),
         pvalue_bin = floor(abs(log_pvalue)*10)/10) %>%
  filter(pvalue_bin >= 3) %>%
  group_by(direction) %>%
  mutate(value = direction == "Positive") %>%
  pull(value) -> outliers.direction
  
prop.test(sum(outliers.direction),length(outliers.direction), 0.5)
dev.off()
 
  
  
all_data_cum %>%
  mutate(window = floor(start/10000)*10000) %>%
  group_by(chr_n,window) %>%
  summarize(high_count = sum(log_pvalue > p_threshold)/n()) %>%
  filter(high_count > 0.05) %>%
  mutate(start = window, end = window+10000) %>%
  select(-window) %>%
  inner_join(chr_names) %>% 
  write_tsv(.,"../../data/gene_copies/sebastes_pgls_gwas_positive_windows.txt")

#Close up of candidates
genes <- read_tsv("../../data/meta/Final_filt_BRK_Sebastes_aleutianus.nameedit.genes.txt",
                  col_names = c("scaf","gene_start","gene_end","ID")) %>%   inner_join(chr_names) 
all_data_cum %>%
  mutate(window = floor(start/10000)*10000) %>%
  group_by(chr_n,window) %>%
  summarize(high_count = sum(log_pvalue > p_threshold)/n()) %>%
  filter(high_count > 0.13) %>%
  mutate(start = window, end = window+10000) %>%
  select(-window) -> top_hits

pdf("plots/read_depth_pgls_gwas.tophits.v1.pdf",height=8,width=8)
for (i in seq(10)){
  print(
  all_data_cum %>%
    filter(chr_n == top_hits$chr_n[i],
           start > top_hits$start[i]-50000,
           start < top_hits$end[i]+50000) %>%
    ggplot() +
    geom_point( aes(x=start, y=abs(log_pvalue)), alpha=0.8, size=1.3,color=colors[2]) +
    theme_cowplot() +
    geom_segment(data=genes %>%
                   filter(chr_n == top_hits$chr_n[i],
                          gene_start > top_hits$start[i]-50000,
                          gene_end < top_hits$end[i]+50000),
                 aes(x=gene_start,xend=gene_end,y=-0.1,yend=-0.1),color="black",size=2) +
    xlab("Bp") + 
    ylab("Log p-value") +
    ggtitle(top_hits$chr_n[i])
  )
}
dev.off()
  
#Saving data
all_data_cum %>%
  select(-tot, -log_pvalue) %>%
  rename(start_cumulative = start_cum) %>%
  write_tsv("../../data/gene_copies/read_depth_pgls_gwas.txt.gz")
  

all_data_cum %>%
  mutate(window = floor(start/10000)*10000) %>%
  group_by(chr_n,window) %>%
  summarize(mean_log_p = median(abs(log_pvalue)),
            high_count = sum(log_pvalue > p_threshold)/n()) -> positive_windows
#Other peak Chr03  6120000-6130000
