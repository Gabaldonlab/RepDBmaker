# Load required libraries
suppressMessages(library(tidyverse))
library(patchwork)

theme_set(theme_classic())

color_db <- c("#c93f55", "#eacc62", "#469d76", "#924099")
names(color_db) <- c("uniprot", "eukprot", "p10k", "custom")

# ---- selection parameters (configurable via dbs.build.repdb.select) ---------
remove_duplicated_species <- as.logical(snakemake@params[["remove_duplicated_species"]])
top_n_genuses <- as.integer(snakemake@params[["top_n_genuses"]])
# clade-reduction table, passed as three parallel vectors
reduce_ranks <- as.character(snakemake@params[["reduce_ranks"]])
reduce_taxa <- as.character(snakemake@params[["reduce_taxa"]])
reduce_ns <- as.integer(snakemake@params[["reduce_ns"]])

# full rank names (or single-letter codes) -> the lineage columns used below
rank_to_col <- c(phylum = "p", class = "c", order = "o", family = "f",
                 p = "p", c = "c", o = "o", f = "f")

# Keep the top-n most complete genomes per family within a clade defined by
# (rank == taxon). Genomes with no resolved family are all kept and are NOT
# subject to the cap: get_tax.R fills unplaced families as "unassigned_<ancestor>"
# (and eukprot/custom may still leave them empty), and we don't know where they
# belong taxonomically, so it makes no sense to rank them within a family.
reduce_clade <- function(data, rank, taxon, n) {
  col <- rank_to_col[[tolower(rank)]]
  if (is.null(col)) {
    stop(paste0("reduce_abundant_clades: unknown rank '", rank, "'"))
  }
  clade <- data[data[[col]] == taxon, ]
  no_family <- clade$f == "" | grepl("^unassigned_", clade$f)
  capped <- clade[!no_family, ] %>%
    group_by(k, p, c, o, f) %>%
    slice_max(completeness, n = n, with_ties = FALSE) %>%
    ungroup()
  bind_rows(capped, clade[no_family, ])
}

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


# first of all remove duplicated species (most complete proteome per species)
if (isTRUE(remove_duplicated_species)) {
    df <- df %>%
        group_by(k, p, c, o, f, g, s) %>%
        slice_max(completeness, n = 1, with_ties = FALSE) %>%
      ungroup()
    paste("After removing duplicated species there are", nrow(df), "genomes")
} else {
    paste("Keeping duplicated species:", nrow(df), "genomes")
}


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
  slice_max(completeness, n = top_n_genuses, with_ties = FALSE) %>%
  ungroup()
# those with no genus are all included
no_genus <- filter(df, g=="")

df <- filter(df, mnemo %in% union(reduced_df$mnemo, no_genus$mnemo))

paste("After keeping max", top_n_genuses, "species per genus there are", nrow(df), "genomes")

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

# reduce over-represented clades: for each configured (rank, taxon, n), keep at
# most n genomes per family within that clade (see reduce_clade() above).
if (length(reduce_ranks) > 0) {
    # rows belonging to at least one configured clade (to drop before adding
    # back the reduced sets); computed on df before any filtering.
    in_clade <- Reduce(`|`, Map(function(rank, taxon) {
        df[[rank_to_col[[tolower(rank)]]]] == taxon
    }, reduce_ranks, reduce_taxa), rep(FALSE, nrow(df)))

    reduced <- bind_rows(Map(function(rank, taxon, n) {
        reduce_clade(df, rank, taxon, n)
    }, reduce_ranks, reduce_taxa, reduce_ns))

    df <- df %>%
        mutate(.in_clade = in_clade) %>%
        filter(!mnemo %in% manual_keep$mnemo, !.in_clade) %>%
        select(-.in_clade) %>%
        bind_rows(reduced, manual_keep, custom_keep) %>%
        distinct(mnemo, .keep_all = TRUE)

    reduce_labels <- paste0(reduce_taxa, " (", reduce_ranks, ", n=", reduce_ns, ")")
    paste("After reducing", paste(reduce_labels, collapse = "; "),
          "and using all the genomes from", paste(to_keep, collapse = ","),
          "there are", nrow(df), "genomes")
} else {
    df <- df %>%
        filter(!mnemo %in% manual_keep$mnemo) %>%
        bind_rows(manual_keep, custom_keep) %>%
        distinct(mnemo, .keep_all = TRUE)
    paste("No clade reduction; there are", nrow(df), "genomes")
}

final <- df %>%
    # filter(db!="p10k" | o != "Ciliophora") %>% 
    group_by(p, c, db) %>% 
    count() %>% 
    ggplot(aes(n, c, fill=db)) +
    labs(subtitle = paste("Reduced clades:", nrow(df))) +
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