library(tidyverse)

files <- list.files("../../data/gene_copies/aleutianus_genome/full_genome_scan_results/", "var")
chr_list <- read_tsv("../../data/meta/chromosome_numbers.txt")
all_variance <- tibble()
for (file in files){
  tmp <- read_tsv(paste0("../../data/gene_copies/aleutianus_genome/full_genome_scan_results/", file))
  all_variance <- rbind(all_variance,tmp)
}
files <- list.files("../../data/gene_copies/aleutianus_genome/full_genome_scan_results/", "74__15")
data <- read_tsv(paste0("../../data/gene_copies/aleutianus_genome/full_genome_scan_results/",files[1])) %>%
  rename(scaffold = scaf)

all_variance %>%
  filter(scaffold == "PGA_scaffold_74__15") %>%
  filter(start > 37600000, start <37800000) %>%
  inner_join(data) %>%
  filter(-log10(pvalue) > 3 ) %>% 
  summarize(mean = mean(sd))
  ggplot(aes(sd)) +
  geom_density()
all_variance %>%
  mutate(type = case_when(sd > 1 ~ "high",
                          TRUE ~ "low")) %>%
  group_by(type) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
  ggplot() +
  geom_density(aes(sd))
