# Build the RepDBmaker UNIVERSE: the single enriched artifact that is the seam
# between the two pipelines.
#
#   Pipeline 1 (sampling / curation) PRODUCES the universe.
#   Pipeline 2 (construction)        CONSUMES it (taxdump, selection, db builds).
#
# One row per available organism (mnemo) across ALL sources, carrying everything
# the downstream pipeline needs so it never has to touch the fragile taxonomy
# harmonization again:
#
#   mnemo          organism id (also the fetch key)
#   source_db      gtdb / virus / uniprot / eukprot / p10k / custom
#   k..s           the 7 harmonized ranks (for the taxdump and for selection)
#   data_type      genome / transcriptome / single-cell / SAG / MAG ...
#   completeness   BUSCO / CheckM2 score where one exists (euks + gtdb), else NA
#   annotated      FALSE for known-unannotated proteomes (dropped by selection)
#
# The universe doubles as the "available proteomes" menu and as the base for the
# per-db provenance table (make_db_meta adds only the build-specific columns).
# It is versioned per release: freezing universe.tsv freezes taxonomy + BUSCO,
# so Pipeline 2 reproduces a release deterministically without re-harmonizing.

suppressMessages(library(tidyverse))

opt_in <- function(name) {
  v <- tryCatch(snakemake@input[[name]], error = function(e) NULL)
  if (is.null(v) || length(v) == 0 || !nzchar(v[[1]]) || !file.exists(v[[1]])) NULL else v[[1]]
}

norm_data_type <- function(x) {
  x <- tolower(replace_na(as.character(x), ""))
  x <- str_replace_all(x, "single[ _]?cell", "single-cell")
  x <- str_replace_all(x, "\\bsag\\b", "single-cell")
  x <- str_replace_all(x, "wgs|whole genome|assembly", "genome")
  na_if(str_squish(x), "")
}

# ---- universe membership + harmonized ranks (id + 7 ranks) -------------------
tax <- read_delim(snakemake@input[["tax"]], delim = "\t",
  col_names = c("mnemo", "k", "p", "c", "o", "f", "g", "s"),
  col_types = cols(.default = "c"), show_col_types = FALSE) %>%
  # id prefix identifies the tracked sources first, so a non-eukaryotic custom
  # (CUS) proteome stays "custom"; the accession-id sources fall back to the
  # superkingdom (which disambiguates GCF: GTDB prokaryote vs NCBI virus).
  mutate(source_db = case_when(
    str_starts(mnemo, "UP")         ~ "uniprot",
    str_starts(mnemo, "EP")         ~ "eukprot",
    str_starts(mnemo, "P10")        ~ "p10k",
    str_starts(mnemo, "CUS")        ~ "custom",
    k %in% c("Bacteria", "Archaea") ~ "gtdb",
    k == "Viruses"                  ~ "virus",
    TRUE                            ~ "other"
  ))

# ---- per-source completeness + data_type + annotation status ----------------
per_source <- list()

per_source$uniprot <- read_delim(snakemake@input[["up_stats"]], show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  transmute(mnemo = proteome_id,
            completeness = parse_number(gsub("\\[.*", "", busco)),
            data_type = "genome", annotated = TRUE)

ep_busco <- read_delim(snakemake@input[["ep_stats"]], show_col_types = FALSE) %>%
  transmute(mnemo = gsub("_.*", "", Input_file), completeness = Complete)
ep_type <- read_delim(snakemake@input[["ep"]], show_col_types = FALSE) %>%
  transmute(mnemo = EukProt_ID, data_type = norm_data_type(Data_Source_Type))
per_source$eukprot <- left_join(ep_busco, ep_type, by = "mnemo") %>% mutate(annotated = TRUE)

# P10K: n_genes == -1 marks unannotated proteomes (selection drops them)
per_source$p10k <- read_delim(snakemake@input[["p10k_stats"]], show_col_types = FALSE) %>%
  transmute(mnemo = p10k_id, completeness = completeness,
            data_type = norm_data_type(strategy), annotated = n_genes != -1)

custom_type <- read_delim(snakemake@input[["custom_table"]], show_col_types = FALSE) %>%
  transmute(mnemo = ID, data_type = norm_data_type(Data_type), annotated = TRUE)
custom_busco_path <- opt_in("custom_busco")
if (!is.null(custom_busco_path)) {
  cb <- read_delim(custom_busco_path, show_col_types = FALSE)
  cb <- tibble(mnemo = cb[[1]], completeness = suppressWarnings(parse_number(as.character(cb[[2]]))))
  per_source$custom <- left_join(custom_type, cb, by = "mnemo")
} else {
  per_source$custom <- mutate(custom_type, completeness = NA_real_)
}

gtdb_meta_path <- opt_in("gtdb_meta")
if (!is.null(gtdb_meta_path)) {
  gm <- read_delim(gtdb_meta_path, show_col_types = FALSE)
  comp_col <- intersect("checkm2_completeness", names(gm))[1]
  type_col <- intersect("ncbi_genome_category", names(gm))[1]
  per_source$gtdb <- tibble(
    mnemo = gsub("_", "", substr(gm$accession, 4, 16)),
    completeness = if (!is.na(comp_col)) suppressWarnings(as.numeric(gm[[comp_col]])) else NA_real_,
    data_type = if (!is.na(type_col)) norm_data_type(gm[[type_col]]) else "genome",
    annotated = TRUE
  )
  per_source$gtdb$data_type[is.na(per_source$gtdb$data_type)] <- "genome"
}

annot <- bind_rows(per_source)

# ---- assemble the universe ---------------------------------------------------
universe <- tax %>%
  left_join(annot, by = "mnemo") %>%
  mutate(annotated = replace_na(annotated, TRUE)) %>%
  # viruses have no completeness concept; keep NA rather than a misleading value
  arrange(mnemo) %>%
  select(mnemo, source_db, k, p, c, o, f, g, s,
         data_type, completeness, annotated)

write_delim(universe, snakemake@output[[1]], delim = "\t")

message(sprintf("universe: %d organisms across %d sources -> %s",
                nrow(universe), n_distinct(universe$source_db), snakemake@output[[1]]))
