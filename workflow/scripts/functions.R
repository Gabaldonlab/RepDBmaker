# Gini coefficient: 0 = every class equally represented, 1 = all sequences in a
# single class.
#
# The previous implementation (adapted from
# https://raw.githubusercontent.com/Arcadia-Science/prehgt/ba092b0e01a20688de56292883d7e35be0e2db08/bin/blastp_to_hgt_candidates_kingdom.R)
# returned 1 - G: it gave 1.000 for a perfectly even vector and 0.020 for an
# extremely skewed one, i.e. an evenness index with the axis label of a Gini
# index. Anything comparing "taxonomic balance" across databases with it was
# inverted. Verify with: gini(rep(1, 10)) == 0.
gini <- function(x) {
  x_sorted <- sort(x)
  n <- length(x_sorted)
  2 * sum(seq_len(n) * x_sorted) / (n * sum(x_sorted)) - (n + 1) / n
}

dmnd_cols <- c(
  "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send",
  "evalue", "bitscore", "qcovhsp", "qlen", "slen"
)

# --- taxonomy-table helpers (used by the harmonization QC report) -------------
# these expect the tidyverse to be attached when they are called.

RANKS <- c("d", "p", "c", "o", "f", "g", "s")

# resolve a path whether run from the repo root or a notebook subdir
rp <- function(p) if (!file.exists(p) && file.exists(file.path("..", "..", p))) file.path("..", "..", p) else p

# read a "mnemo <tab> d__..;..;s__.." table into mnemo + 7 rank cols (+ raw
# lineage), dropping identical duplicate rows
read_prefixed <- function(path) {
  read_tsv(rp(path), col_names = c("mnemo", "lineage"),
           col_types = cols(.default = "c"), name_repair = "minimal") %>%
    distinct() %>%
    separate(lineage, RANKS, sep = ";", fill = "right", extra = "drop", remove = FALSE) %>%
    mutate(across(all_of(RANKS), ~ trimws(gsub("^[a-z]__", "", .))))
}

# TRUE where a and b are both non-empty and differ
changed <- function(a, b) !is.na(a) & !is.na(b) & a != "" & b != "" & a != b

# mnemos that carry more than one distinct lineage
ambiguous_mnemos <- function(harm) harm %>% distinct(mnemo, lineage) %>%
  count(mnemo) %>% filter(n > 1) %>% pull(mnemo)

# for each `child` value, keep its single most frequent `parent` (used to force
# a union of taxonomies into a valid, single-parent hierarchy)
modal_parent <- function(df, child, parent) df %>%
  count(.data[[child]], .data[[parent]], name = "nn") %>%
  group_by(.data[[child]]) %>% slice_max(nn, n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(all_of(c(child, parent)))