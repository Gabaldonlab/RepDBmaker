library(tidyverse)

tax <- read_delim(snakemake@input[["tax"]], 
           col_names = c("file", "k", "p", "c", "o", "f", "g", "s"))
stats <- read_delim(snakemake@input[["stats"]]) %>% 
  left_join(tax) %>% 
  rename("mnemo"="file")

up_stats <- read_delim(snakemake@input[["up_stats"]], show_col_types = F) %>% 
  janitor::clean_names() %>% 
  mutate(completeness=parse_number(gsub("\\[.*", "", busco))) %>% 
  select(proteome_id, completeness) %>% 
  mutate(data_type="genome")
colnames(up_stats) <- c("mnemo", "completeness", "data_type")
ep_stats <- read_delim(snakemake@input[["ep_stats"]], show_col_types = F) %>% 
  mutate(mnemo=gsub("_.*", "", Input_file)) %>% 
  select(mnemo, Complete)
colnames(ep_stats) <- c("mnemo", "completeness")
ep_source <- read_delim(snakemake@input[["ep"]]) %>% 
  select(EukProt_ID, Data_Source_Type)
colnames(ep_source) <- c("mnemo", "data_type")
ep_stats <- left_join(ep_stats, ep_source)

p10k_stats <- read_delim(snakemake@input[["p10k_stats"]], show_col_types = F) %>% 
	filter(n_genes!=-1) %>%
	select(p10k_id, completeness, strategy)
colnames(p10k_stats) <- c("mnemo", "completeness", "data_type")

custom_stats <- read_delim(snakemake@input[["custom_busco"]]) %>% 
  select(file, complete)
colnames(custom_stats) <- c("mnemo", "completeness")
custom_source <- read_delim(snakemake@input[["custom_table"]]) %>% 
  select(ID, Data_type)
colnames(custom_source) <- c("mnemo", "data_type")
custom_stats <- left_join(custom_stats, custom_source)

contaminants <- table(str_split(readLines(snakemake@input[["contaminants"]]), "_", simplify = T)[,2]) %>% 
  data.frame()
colnames(contaminants) <- c("mnemo", "n_contaminants")

final_stats <- stats %>% 
  filter(k=="Eukaryota") %>% 
  left_join(rbind(up_stats, ep_stats, p10k_stats, custom_stats)) %>% 
  left_join(contaminants) %>% 
  mutate(data_type=gsub("Transcriptome", "transcriptome", 
                        gsub("Genome", "genome",
                             gsub("Single cell", "single-cell", data_type)))) %>% 
  mutate(prop_contaminants = n_contaminants/num_seqs) %>% 
  select(mnemo, k, p, c, o, f, g, s, data_type, completeness, contains("contaminants"),
         num_seqs, sum_len, min_len, avg_len, max_len)

write_delim(final_stats, snakemake@output[[1]], delim = "\t")

