#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))

gtdb <- read_delim(stdin(), show_col_types = FALSE)

gtdb %>% 
    filter(gtdb_representative) %>% 
    # filter(checkm_completeness>50) %>% 
    separate(gtdb_taxonomy, c("k", "p", "c", "o", "f", "g", "s"), ";") %>% 
    group_by(k,p,c,o,f,g) %>% 
    slice_max(checkm_completeness, with_ties = FALSE) %>% 
    write_delim(stdout(), delim = "\t")
