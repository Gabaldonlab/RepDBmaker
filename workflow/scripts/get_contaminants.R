suppressPackageStartupMessages(library(tidyverse))

min_prop <- snakemake@params[["prop_euka"]]
# max_size <- snakemake@params[["size_cluster"]]

clusters <- read_delim(snakemake@input[[1]], 
                       col_names = c("rep", "seq", "taxid", "tax"), delim = "\t")

euks <- clusters %>% 
  filter(grepl("Eukaryota;", tax)) %>% 
  separate(tax, c("k", "p", "c", "f", "o", "g", "s"), ";")

single_euk <- euks %>% 
  group_by(rep) %>% 
  count() %>% 
  filter(n==1)

clusters_nonsingle <- filter(clusters, !rep %in% single_euk$rep)
rm(clusters)

df_euka_prop <- clusters_nonsingle %>% 
  separate(tax, c("k", "p", "c", "f", "o", "g", "s"), ";") %>% 
  group_by(rep) %>%
  mutate(n_clu=n(), n_sps=n_distinct(s)) %>% 
  group_by(rep, k, n_clu, n_sps) %>% 
  count() %>% 
  pivot_wider(names_from = k, values_from = n) %>% 
  mutate(euka_prop = Eukaryota/n_clu) %>% 
  ungroup()

clusters_to_filter <- filter(df_euka_prop, euka_prop <= min_prop)$rep #  | n_clu<max_size

euks %>%
  filter(rep %in% c(clusters_to_filter, single_euk$rep)) %>% 
  pull(seq) %>% 
  writeLines(snakemake@output[[1]])
