suppressPackageStartupMessages(library(tidyverse))

####### Prokaryotes ######
gtdb_tax <- read_delim(snakemake@input[["gtdb_tax"]],
  col_names = c("mnemo", "lineage"),
  show_col_types = FALSE, delim = "\t"
)

gtdb <- read_delim(snakemake@input[["gtdb"]],
  show_col_types = FALSE, delim = "\t"
) %>%
  select(
    accession,
    checkm2_completeness,
    # checkm2_contamination,
    ncbi_assembly_level
  ) %>%
  mutate(mnemo = gsub("^GB|^RS|\\..*", "", gsub("_", "", accession))) %>%
  rename(
    "Data_type" = "ncbi_assembly_level",
    "Completeness" = "checkm2_completeness"
  ) %>%
  select(mnemo, Completeness, Data_type) %>%
  left_join(gtdb_tax)

######### Viruses ########
virus_tax <- read_delim(snakemake@input[["viruses_tax"]],
  col_names = c("mnemo", "lineage"),
  show_col_types = FALSE, delim = "\t"
)

virus_meta <- read_delim(snakemake@input[["viruses"]],
  col_names = FALSE, show_col_types = FALSE, delim = "\t"
) %>%
  select(X1, X12) %>%
  mutate(
    mnemo = gsub("\\..*", "", gsub("_", "", X1)),
    Completeness = NA
  ) %>%
  rename("Data_type" = "X12") %>%
  select(mnemo, Data_type, Completeness) %>%
  # here you need to right join as 100s viruses are in NCBI viruses
  # but return error in nuccore (GCF_000857205 for example)
  right_join(virus_tax)

######### Eukaryotes ######

euka_tax <- read_delim(snakemake@input[["euka_tax"]],
  col_names = c("mnemo", "lineage"),
  show_col_types = FALSE, delim = "\t"
)
# EP
ep <- read_delim(snakemake@input[["ep_included"]], show_col_types = F)
ep_stats <- read_delim(snakemake@input[["ep_busco"]], show_col_types = F) %>%
  mutate(mnemo = gsub("_.*", "", Input_file)) %>%
  select(mnemo, Complete) %>%
  left_join(select(ep, c(EukProt_ID, Data_Source_Type)),
    by = c("mnemo" = "EukProt_ID")
  )
colnames(ep_stats) <- c("mnemo", "Completeness", "Data_type")

p10k_stats <- read_delim(snakemake@input[["p10k"]], show_col_types = F)
# some p10k genomes are not annotated! Remove them
p10k_unannotated <- p10k_stats[p10k_stats$n_genes == -1, ]$p10k_id
p10k_stats_red <- select(p10k_stats, p10k_id, completeness, strategy) %>%
  mutate(strategy = gsub(
    "single cell", "single-cell",
    gsub("\\ \\(.*", "", str_to_lower(strategy))
  )) %>%
  filter(!p10k_id %in% p10k_unannotated)
colnames(p10k_stats_red) <- c("mnemo", "Completeness", "Data_type")

up_stats <- read_delim(snakemake@input[["up"]], show_col_types = F) %>%
  janitor::clean_names() %>%
  mutate(completeness = parse_number(gsub("\\[.*", "", busco))) %>%
  select(proteome_id, completeness) %>%
  mutate(Data_type = "genome")
colnames(up_stats) <- c("mnemo", "Completeness", "Data_type")

# On hold as the user should now what they have available
# cus_stats <- read_delim("../../DB/new_genomes/results/stats/custom_busco.tsv") %>%
#   mutate(strategy = NA) %>%
#   select(file, complete, strategy)
# colnames(cus_stats) <- c("mnemo", "Completeness", "Strategy")

df <- bind_rows(up_stats, ep_stats, p10k_stats_red) %>%
  left_join(euka_tax) %>%
  bind_rows(gtdb) %>%
  bind_rows(virus_meta) %>%
  separate(lineage, c("k", "p", "c", "o", "f", "g", "s"), ";") %>%
  mutate_all(~ gsub("[a-z]__", "", .)) %>%
  arrange(k, p, c, o, f, g, s) %>%
  mutate(Data_type = gsub(
    "Est", "EST",
    gsub(
      ",est", "",
      gsub(
        "Single-cell", "S.C.",
        str_to_sentence(Data_type)
      )
    )
  )) %>%
  #  few uniprot genomes fail for taxonomy
  filter(!is.na(k))

write_delim(df, snakemake@output[[1]], delim = "\t")
