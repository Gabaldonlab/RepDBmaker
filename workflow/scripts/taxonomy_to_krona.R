#!/usr/bin/env Rscript

# Convert a RepDBmaker taxonomy table into the text format expected by
# KronaTools' ktImportText: one line per proteome, a magnitude (1) followed by
# the taxonomic ranks from domain to species. Empty ranks become "unassigned" so
# the hierarchy keeps its depth.
#
# Two input formats (param `fmt`):
#   "prefixed"    mnemo <tab> d__..;p__..;..;s__..           (the default; also
#                 the UniProt NCBI taxonomy)
#   "p10k_native" mnemo <tab> species <tab> Super;P_..;C_..  (P10K's native
#                 referable lineage)
# This lets a source' original schema and its UniEuk-harmonized version both be
# charted, so the harmonization can be compared.
#
# Usage (standalone):
#   Rscript taxonomy_to_krona.R <taxonomy.tsv> <out.krona.txt> [fmt]
# Also callable from Snakemake (input[[1]] / output[[1]] / params$fmt).

suppressPackageStartupMessages(library(tidyverse))
source("workflow/scripts/functions.R")   # read_prefixed, RANKS, rp

# P10K native lineage: 3 cols, ranks tagged Super;P_..;C_..;O_..;F_..;G_..;S_..
read_p10k_native <- function(path) {
  read_tsv(rp(path), col_names = c("mnemo", "species", "lineage"),
           col_types = cols(.default = "c"), name_repair = "minimal") %>%
    distinct() %>%
    separate(lineage, RANKS, sep = ";", fill = "right", extra = "drop", remove = FALSE) %>%
    mutate(across(all_of(RANKS), ~ trimws(gsub("^[A-Z]_", "", .))))
}

if (exists("snakemake")) {
  infile  <- snakemake@input[["tax"]]
  outfile <- snakemake@output[[1]]
  fmt     <- snakemake@params[["fmt"]]
  root    <- snakemake@params[["root"]]
} else {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) stop("usage: taxonomy_to_krona.R <taxonomy.tsv> <out.krona.txt> [fmt] [root]")
  infile <- args[[1]]; outfile <- args[[2]]
  fmt  <- if (length(args) >= 3) args[[3]] else "prefixed"
  root <- if (length(args) >= 4) args[[4]] else ""
}
if (is.null(fmt) || length(fmt) == 0 || is.na(fmt) || fmt == "") fmt <- "prefixed"
if (is.null(root) || length(root) == 0 || is.na(root)) root <- ""

tax <- switch(fmt,
  prefixed    = read_prefixed(infile),
  p10k_native = read_p10k_native(infile),
  stop("unknown fmt: ", fmt))

# fill an empty domain with `root` (e.g. the UniProt NCBI table has no
# superkingdom, so root it at "Eukaryota" for a clean comparison)
if (root != "") tax <- tax %>% mutate(d = ifelse(is.na(d) | d == "", root, d))
if (nrow(tax) == 0) {
  # keep an (empty) file so the plot rule still has an input
  file.create(outfile)
} else {
  tax %>%
    mutate(across(all_of(RANKS), ~ ifelse(is.na(.) | . == "", "unassigned", .))) %>%
    transmute(magnitude = 1L, d, p, c, o, f, g, s) %>%
    write_tsv(outfile, col_names = FALSE)
}
