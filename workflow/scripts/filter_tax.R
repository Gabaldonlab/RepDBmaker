# Load required libraries
suppressMessages(library(tidyverse))
library(patchwork)

theme_set(theme_classic())

color_db <- c("#c93f55", "#eacc62", "#469d76", "#924099")
names(color_db) <- c("uniprot", "eukprot", "p10k", "custom")

exclude <- readLines(snakemake@input[["exclude"]])

up_stats <- read_delim(snakemake@input[["up_stats"]], show_col_types = F) %>% 
    janitor::clean_names() %>% 
    mutate(completeness=parse_number(gsub("\\[.*", "", busco))) %>% 
    select(proteome_id, completeness)
colnames(up_stats) <- c("mnemo", "completeness")
ep_stats <- read_delim(snakemake@input[["ep_stats"]], show_col_types = F) %>% 
    mutate(mnemo=gsub("_.*", "", Input_file)) %>% 
    select(mnemo, Complete)
colnames(ep_stats) <- c("mnemo", "completeness")
p10k_stats <- read_delim(snakemake@input[["p10k_stats"]], show_col_types = F)

# some p10k genomes are not annotated! Remove them
p10k_unannotated <- p10k_stats[p10k_stats$n_genes==-1, ]$p10k_id

paste("There are ", length(p10k_unannotated), "P10K unannotated genomes")

p10k_stats_red <- select(p10k_stats, p10k_id, completeness)
colnames(p10k_stats_red) <- c("mnemo", "completeness")


og <- read_delim(c(snakemake@input[["up"]], snakemake@input[["ep"]], 
                   snakemake@input[["p10k"]], snakemake@input[["custom"]]),
                 col_names = c("mnemo", "lineage"), delim = "\t", show_col_types = F)
df <- og %>% 
    separate(lineage, c("k", "p", "c", "o", "f", "g", "s"), ";") %>% 
    mutate_all(~gsub("[a-z]__", "", .)) %>% 
    rowwise() %>% 
    mutate(db = case_when(grepl("^UP", mnemo) ~ "uniprot",
                          grepl("^EP", mnemo) ~ "eukprot",
                          grepl("^P10", mnemo) ~ "p10k",
                          grepl("^CUS", mnemo) ~ "custom"),
           db = factor(db, levels = c("custom","uniprot", "eukprot", "p10k")))  %>% 
    left_join(rbind(up_stats, ep_stats, p10k_stats_red)) %>% 
    filter(!mnemo %in% p10k_unannotated) %>% 
    filter(!mnemo %in% exclude)



# keep all custom genomes
custom_keep <- df %>%
  filter(db=="custom")
paste("There are", nrow(custom_keep), "custom genomes")

df <- df %>%
  filter(db!="custom")

paste("There are", nrow(df), "genomes")
# keep all:
# - Crums
# - Telonemia
# - Ancyro
# - Malawi
# - Anaeromaebae??
to_keep <- readLines(snakemake@input[["to_keep"]])
manual_keep <- df %>%
    filter(rowSums(across(everything(), ~. %in% to_keep)) > 0)
    # filter(p=="CRuMs" | p=="Ancyromonadida" | p=="Malawimonadida" | 
    #        c=="Telonemia" | c=="Rhodelphidia" | c=="Anaeramoebidae")
paste(nrow(manual_keep), "genomes were forcefully kept")

raw <- df %>%
    group_by(p, c, db) %>% 
    count() %>% 
    ggplot(aes(n, c, fill=db)) +
    labs(subtitle = paste("Raw:", nrow(df))) +
    facet_grid(p~., scales = "free", space = "free") +
    geom_bar(stat = "identity", position = "dodge") + 
    geom_text(aes(label=n), hjust=0, position = position_dodge(width = .9), size=2) +
    scale_fill_manual(values=color_db)


# first of all remove duplicated species
df <- df %>% 
    group_by(k, p, c, o, f, g, s) %>% 
    arrange(desc(completeness), db) %>% 
    slice(1)

paste("After removing duplicated species there are", nrow(df), "genomes")


no_dup_sps <- df %>%
    group_by(p, c, db) %>% 
    count() %>% 
    ggplot(aes(n, c, fill=db)) +
    labs(subtitle = paste("No dup sp:", nrow(df))) +
    facet_grid(p~., scales = "free", space = "free") +
    geom_bar(stat = "identity", position = "dodge") + 
    geom_text(aes(label=n), hjust=0, position = position_dodge(width = .9), size=2) +
    scale_fill_manual(values=color_db) +
    theme(axis.text.y = element_blank())

# genuses that have more than one representative get only the most complete
reduced_df <- df %>% 
    group_by(k, p, c, o, f, g) %>% 
    arrange(desc(completeness), db) %>% 
    slice(1)
# those with no genus are all included
no_genus <- filter(df, g=="")

df <- filter(df, mnemo %in% union(reduced_df$mnemo, no_genus$mnemo))

paste("After keeping one species per genus there are", nrow(df), "genomes")


no_dup_genus <- df %>%
    # filter(db!="p10k" | o != "Ciliophora") %>% 
    group_by(p, c, db) %>% 
    count() %>% 
    ggplot(aes(n, c, fill=db)) +
    labs(subtitle = paste("1 x genus:", nrow(df))) +
    facet_grid(p~., scales = "free", space = "free") +
    geom_bar(stat = "identity", position = "dodge") + 
    geom_text(aes(label=n), hjust=0, position = position_dodge(width = .9), size=2) +
    scale_fill_manual(values=color_db) +
    theme(axis.text.y = element_blank())

opi <- df %>% 
    filter(c=="Opisthokonta") %>% 
    group_by(k, p, c, o, f) %>% 
    arrange(desc(completeness), db) %>% 
    slice_head(n = 20)
opi_nofamily <- filter(df, c=="Opisthokonta", f=="") %>% 
    filter(!mnemo %in% opi$mnemo)

ciliates <- df %>% 
    filter(o=="Ciliophora") %>% 
    group_by(k, p, c, o, f) %>% 
    arrange(desc(completeness), db) %>% 
    slice_head(n = 20)
ciliates_nofamily <- filter(df, o=="Ciliophora", f=="") %>% 
    filter(!mnemo %in% ciliates$mnemo)


df <- df %>% 
    filter(!mnemo %in% manual_keep$mnemo) %>% 
    filter(c!="Opisthokonta", o!="Ciliophora") %>% 
    rbind(opi, opi_nofamily) %>% 
    rbind(ciliates, ciliates_nofamily) %>% 
    rbind(manual_keep) %>% 
    rbind(custom_keep)

paste("After keeping 20 genomes per family of ciliates and opisthokonta and using all the genomes from", 
      paste(to_keep, collapse=","), "there are", nrow(df), "genomes")

final <- df %>%
    # filter(db!="p10k" | o != "Ciliophora") %>% 
    group_by(p, c, db) %>% 
    count() %>% 
    ggplot(aes(n, c, fill=db)) +
    labs(subtitle = paste("20 x opi and ciliates family:", nrow(df))) +
    facet_grid(p~., scales = "free", space = "free") +
    geom_bar(stat = "identity", position = "dodge") + 
    geom_text(aes(label=n), hjust=0, position = position_dodge(width = .9), size=2) +
    scale_fill_manual(values=color_db) +
    theme(axis.text.y = element_blank())


final_plot <- (raw | no_dup_sps | no_dup_genus | final) + plot_layout(guides = 'collect')

ggsave(snakemake@output[["plot"]], final_plot, width = 12, height = 12)

filter(og, mnemo %in% df$mnemo) %>% 
    write_delim(snakemake@output[["tax"]], delim = "\t", col_names = F)


# maybe also do not remove if they dont have a complete enough proteome

# if there are duplicated proteomes its ok as you are going to cluster for very high identity 
# or you will look for duplicates with seqkit