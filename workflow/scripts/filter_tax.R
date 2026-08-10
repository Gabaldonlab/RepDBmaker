# Load required libraries
suppressMessages(library(tidyverse))
library(patchwork)

theme_set(theme_classic())

color_db <- c("#c93f55", "#eacc62", "#469d76", "#924099")
names(color_db) <- c("uniprot", "eukprot", "p10k", "custom")

# ---- selection parameters (configurable via dbs.build.repdb.select) ---------
remove_duplicated_species <- as.logical(snakemake@params[["remove_duplicated_species"]])
if (is.na(remove_duplicated_species)) {
    print("remove_duplicated_species parameter is not a logical value in config.yaml")
    quit(status = 1)
}
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

# Eukaryote selection operates on the UNIVERSE (id, source, ranks, completeness,
# annotated). Only eukaryotes are downsampled here; prokaryotes/viruses are taken
# wholesale into RepDB later. Reading the universe makes the selection a pure,
# deterministic function of (universe, config) - so it reproduces exactly from a
# pinned universe without re-running the taxonomy harmonization.
EUK_SOURCES <- c("custom", "uniprot", "eukprot", "p10k")
df <- read_delim(snakemake@input[["universe"]], delim = "\t", show_col_types = FALSE,
                 col_types = cols(.default = "c")) %>%
    filter(source_db %in% EUK_SOURCES) %>%
    mutate(db = factor(source_db, levels = EUK_SOURCES),
           completeness = suppressWarnings(as.numeric(completeness)),
           annotated = as.logical(annotated)) %>%
    filter(is.na(annotated) | annotated) %>%   # drop known-unannotated proteomes (P10K n_genes == -1)
    filter(!mnemo %in% exclude) %>%
    # deterministic tie-break: slice_max(with_ties = FALSE) keeps the first row
    # among equal-completeness ties, so fix the order by mnemo (otherwise the
    # selection would depend on the universe row order).
    arrange(mnemo)



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
    scale_fill_manual(values=color_db, drop = FALSE)


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
    scale_fill_manual(values=color_db, drop = FALSE) +
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
    labs(subtitle = paste(top_n_genuses," x genus:", nrow(df))) +
    facet_grid(p~., scales = "free", space = "free") +
    geom_bar(stat = "identity", position = "dodge") + 
    geom_text(aes(label=n), hjust=0, position = position_dodge(width = .9), size=2) +
    scale_fill_manual(values=color_db, drop = FALSE) +
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
    scale_fill_manual(values=color_db, drop = FALSE) +
    theme(axis.text.y = element_blank())


final_plot <- (raw | no_dup_sps | no_dup_genus | final) + plot_layout(guides = 'collect')

ggsave(snakemake@output[["plot"]], final_plot, width = 12, height = 12)

# write the selected eukaryotes as mnemo<TAB>prefixed-lineage (repdb_taxonomy
# only needs the id column, but keep the lineage for readability / QC).
b <- function(x) coalesce(x, "")
df %>%
    transmute(mnemo,
              completeness,
              lineage = paste0("d__", b(k), ";p__", b(p), ";c__", b(c), ";o__", b(o),
                               ";f__", b(f), ";g__", b(g), ";s__", b(s))) %>%
    write_delim(snakemake@output[["tax"]], delim = "\t", col_names = FALSE)


# maybe also do not remove if they dont have a complete enough proteome

# if there are duplicated proteomes its ok as you are going to cluster for very high identity 
# or you will look for duplicates with seqkit
