library(tidyverse)

tax <- read_delim("results/taxonomies/repdb.tsv", 
           col_names = c("file", "k", "p", "c", "o", "f", "g", "s"))
stats <- read_delim("results/meta/repdb_stats.tsv") %>% 
  left_join(tax)

