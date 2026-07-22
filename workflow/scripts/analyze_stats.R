# Assemble a per-organism provenance/metadata table for a database (RepDB or any
# custom database), written to results/meta/<db>_meta.tsv.
#
# Reviewer 1, methodological point #3 (source heterogeneity): integrating sources
# with different foundations produces an output that is more heterogeneous than a
# single flat FASTA suggests. Rather than discarding the intermediate state, this
# script propagates it so users can subset the database along the axes that
# matter (source, data type, completeness metric, upstream/RepDBmaker filtering).
#
# One row per organism (mnemo) across ALL sources (not just eukaryotes) with:
#   taxid                the NCBI-style taxid assigned in the database taxdump
#   source_db            gtdb / virus / uniprot / eukprot / p10k / custom
#   k..s                 the (harmonized) 7-rank lineage as stored in the database
#   data_type            genome / transcriptome / single-cell / SAG / MAG ...
#   completeness         a numeric score WHERE ONE EXISTS ...
#   n_contaminants,      RepDBmaker's own decontamination flagging for this
#   prop_contaminants    organism (0 / NA when decontamination was not run)
#   num_seqs..max_len    seqkit size statistics
#
# Optional inputs (custom_busco, gtdb_meta, contaminants) degrade gracefully.

suppressMessages(library(tidyverse))

opt_in <- function(name) {
  v <- tryCatch(snakemake@input[[name]], error = function(e) NULL)
  if (is.null(v) || length(v) == 0 || !nzchar(v[[1]]) || !file.exists(v[[1]])) NULL else v[[1]]
}
strip_ext <- function(x) str_remove(basename(x), "\\.(faa|fasta|fa)(\\.gz)?$")

# ---- taxonomy of the RepDB composition (ID + 7 unprefixed ranks) -------------
tax <- read_delim(snakemake@input[["tax"]], delim = "\t",
  col_names = c("mnemo", "k", "p", "c", "o", "f", "g", "s"),
  col_types = cols(.default = "c"), show_col_types = FALSE) %>%
  # source database: superkingdom disambiguates GCF (GTDB prok vs NCBI virus)
  mutate(source_db = case_when(
    k %in% c("Bacteria", "Archaea") ~ "gtdb",
    k == "Viruses"                  ~ "virus",
    str_starts(mnemo, "UP")         ~ "uniprot",
    str_starts(mnemo, "EP")         ~ "eukprot",
    str_starts(mnemo, "P10")        ~ "p10k",
    str_starts(mnemo, "CUS")        ~ "custom",
    TRUE                            ~ "other"
  ))

# ---- taxid assigned in the taxdump (taxonkit create-taxdump -A1 taxid.map) ----
# taxid.map maps each accession (= mnemo, the all_taxonomy first column) to its
# taxid, so this is the exact taxid the DIAMOND/MMseqs/BLAST indices use.
taxdump_dir <- opt_in("taxdump")
taxid_map <- if (!is.null(taxdump_dir)) file.path(taxdump_dir, "taxid.map") else ""
if (nzchar(taxid_map) && file.exists(taxid_map)) {
  taxids <- read_delim(taxid_map, delim = "\t", col_names = c("mnemo", "taxid"),
    col_types = cols(.default = "c"), show_col_types = FALSE) %>%
    distinct(mnemo, .keep_all = TRUE)
} else {
  message("taxid.map not found in the taxdump directory; taxid will be NA")
  taxids <- tibble(mnemo = character(0), taxid = character(0))
}


norm_data_type <- function(x) {
  x <- tolower(replace_na(as.character(x), ""))
  x <- str_replace_all(x, "single[ _]?cell", "single-cell")
  x <- str_replace_all(x, "\\bsag\\b", "single-cell")
  x <- str_replace_all(x, "wgs|whole genome|assembly", "genome")
  na_if(str_squish(x), "")
}

# ---- per-source completeness + data_type ------------------------------------
per_source <- list()

# UniProt: reference proteomes (BUSCO), all genome-derived
per_source$uniprot <- read_delim(snakemake@input[["up_stats"]], show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  transmute(mnemo = proteome_id,
            completeness = parse_number(gsub("\\[.*", "", busco)),
            data_type = "genome")

# EukProt: BUSCO score + declared source type
ep_busco <- read_delim(snakemake@input[["ep_stats"]], show_col_types = FALSE) %>%
  transmute(mnemo = gsub("_.*", "", Input_file), completeness = Complete)
ep_type <- read_delim(snakemake@input[["ep"]], show_col_types = FALSE) %>%
  transmute(mnemo = EukProt_ID, data_type = norm_data_type(Data_Source_Type))
per_source$eukprot <- left_join(ep_busco, ep_type, by = "mnemo")

# P10K: completeness + sequencing strategy as data type
per_source$p10k <- read_delim(snakemake@input[["p10k_stats"]], show_col_types = FALSE) %>%
  filter(n_genes != -1) %>%
  transmute(mnemo = p10k_id, completeness = completeness,
            data_type = norm_data_type(strategy))

# custom: declared data type (+ optional BUSCO if the user supplied one)
custom_type <- read_delim(snakemake@input[["custom_table"]], show_col_types = FALSE) %>%
  transmute(mnemo = ID, data_type = norm_data_type(Data_type))
custom_busco_path <- opt_in("custom_busco")
if (!is.null(custom_busco_path)) {
  cb <- read_delim(custom_busco_path, show_col_types = FALSE)
  cb <- tibble(mnemo = cb[[1]], completeness = suppressWarnings(parse_number(as.character(cb[[2]]))))
  per_source$custom <- left_join(custom_type, cb, by = "mnemo")
} else {
  per_source$custom <- mutate(custom_type, completeness = NA_real_)
}

# GTDB: CheckM2 completeness + genome category (optional; accession -> RepDB id)
gtdb_meta_path <- opt_in("gtdb_meta")
if (!is.null(gtdb_meta_path)) {
  gm <- read_delim(gtdb_meta_path, show_col_types = FALSE)
  comp_col <- intersect("checkm2_completeness", names(gm))[1]
  type_col <- intersect("ncbi_assembly_level", names(gm))[1]
  per_source$gtdb <- tibble(
    mnemo = gsub("_", "", substr(gm$accession, 4, 16)),
    completeness = if (!is.na(comp_col)) suppressWarnings(as.numeric(gm[[comp_col]])) else NA_real_,
    data_type = if (!is.na(type_col)) norm_data_type(gm[[type_col]]) else "genome"
  )
  per_source$gtdb$data_type[is.na(per_source$gtdb$data_type)] <- "genome"
}

annot <- bind_rows(per_source)  # mnemo -> completeness, data_type

# ---- RepDBmaker decontamination flagging (optional) --------------------------
# contaminant sequence IDs -> organism = the token before the first "_"
# (same convention as make_db_map / get_mixed_clusters).
contam_path <- opt_in("contaminants")
if (!is.null(contam_path)) {
  ids <- readLines(contam_path)
  ids <- ids[nzchar(ids)]
  contaminants <- tibble(mnemo = sub("_.*", "", ids)) %>%
    count(mnemo, name = "n_contaminants")
} else {
  contaminants <- tibble(mnemo = character(0), n_contaminants = integer(0))
}

# ---- seqkit size statistics for every organism ------------------------------
stats <- read_delim(snakemake@input[["stats"]], show_col_types = FALSE)
names(stats)[1] <- "file"
stats <- stats %>% mutate(mnemo = strip_ext(file)) %>% select(-file)

# ---- assemble ----------------------------------------------------------------
meta <- tax %>%
  left_join(taxids, by = "mnemo") %>%
  left_join(annot, by = "mnemo") %>%
  left_join(stats, by = "mnemo") %>%
  left_join(contaminants, by = "mnemo") %>%
  mutate(
    n_contaminants = replace_na(n_contaminants, 0L),
    prop_contaminants = if_else(!is.na(num_seqs) & num_seqs > 0,
                                n_contaminants / num_seqs, NA_real_),
    # sources with no completeness concept keep NA rather than a misleading 0
    completeness = if_else(completeness_metric == "none", NA_real_, completeness)
  ) %>%
  select(mnemo, taxid, source_db,
         k, p, c, o, f, g, s,
         data_type, completeness,
         n_contaminants, prop_contaminants,
         any_of(c("num_seqs", "sum_len", "min_len", "avg_len", "max_len")))

write_delim(meta, snakemake@output[[1]], delim = "\t")

message(sprintf("%s: %d organisms across %d sources written",
                snakemake@output[[1]], nrow(meta), n_distinct(meta$source_db)))
