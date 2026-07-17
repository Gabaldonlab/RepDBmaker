suppressPackageStartupMessages(library(tidyverse))

# Identify eukaryotic sequences that look like contamination from mixed clusters
# (a cluster containing eukaryotic and non-eukaryotic members). Two flag paths:
#   single_euk    - a lone eukaryote in a mixed cluster (threshold-independent)
#   low_euka_prop - eukaryotes are a minority of the cluster (<= prop_euka)
# Writes, regardless of the downstream hard/soft filter mode:
#   - a data frame of every flagged protein with its cluster properties (`df`)
#   - the plain list of flagged protein IDs (`ids`), used by the filter step

min_prop <- snakemake@params[["prop_euka"]]

ranks <- c("k", "p", "c", "o", "f", "g", "s")

clusters <- read_delim(snakemake@input[[1]],
  col_names = c("rep", "seq", "taxid", "tax"), delim = "\t",
  show_col_types = FALSE
) %>%
  separate(tax, ranks, sep = ";", fill = "right", remove = FALSE)

# per-cluster properties
cluster_props <- clusters %>%
  group_by(rep) %>%
  summarise(
    n_clu = n(),
    n_euka = sum(k == "Eukaryota", na.rm = TRUE),
    n_species = n_distinct(s),
    euka_prop = n_euka / n_clu,
    .groups = "drop"
  )

single_reps <- cluster_props %>% filter(n_euka == 1) %>% pull(rep)
low_prop_reps <- cluster_props %>% filter(n_euka > 1, euka_prop <= min_prop) %>% pull(rep)

flagged <- bind_rows(
  tibble(rep = single_reps, flag = "single_euk"),
  tibble(rep = low_prop_reps, flag = "low_euka_prop")
)

# one row per flagged (eukaryotic) protein, with the cluster context
contaminants <- clusters %>%
  filter(k == "Eukaryota", rep %in% flagged$rep) %>%
  left_join(flagged, by = "rep") %>%
  left_join(cluster_props, by = "rep") %>%
  select(seq, rep, flag, n_clu, n_euka, n_species, euka_prop, all_of(ranks)) %>%
  arrange(rep, seq)

write_tsv(contaminants, snakemake@output[["df"]])
writeLines(contaminants$seq, snakemake@output[["ids"]])
