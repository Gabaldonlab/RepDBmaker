#!/usr/bin/env Rscript

# Validate a RepDBmaker custom-proteome table before running the pipeline.
#
# Custom proteomes are added through the tab-separated table pointed to by
# `files.new_genomes` in the config (e.g. resources/custom_genomes_repdb.csv).
# The pipeline makes several *implicit* assumptions about that table; when they
# are violated the offending rows are not rejected but silently dropped or
# mislabelled downstream (reviewer methodological concern #1). This script makes
# those assumptions explicit and checks them up front.
#
# Required columns (tab-separated, with a header):
#   ID        unique mnemonic, must start with 'CUS'
#   Species   organism name (single organism per row)
#   Data_type e.g. genome / transcriptome (informational)
#   Fasta     path to the proteome FASTA for this entry
#   Lineage   exactly 7 ';'-separated ranks with the prefixes
#             d__ p__ c__ o__ f__ g__ s__ , e.g.
#             d__Eukaryota;p__Discoba;c__Jakobida;o__Jakobida;f__Ophirinina;g__Agogonia;s__Agogonia voluta
#
# Checks (ERROR fails validation / non-zero exit; WARNING is reported only):
#   ERROR   missing required column
#   ERROR   ID does not start with 'CUS'   (else silently dropped by filter_tax.R)
#   ERROR   duplicate ID
#   ERROR   lineage not exactly 7 ranks    (else separate() mangles it)
#   ERROR   lineage rank has wrong/missing prefix (d__,p__,c__,o__,f__,g__,s__)
#   ERROR   lineage rank empty after its prefix
#   ERROR   empty Fasta path
#   ERROR   duplicate Fasta path
#   ERROR   Fasta path missing or empty     (only with --check-fasta)
#   ERROR   domain (d__) is not a recognized superkingdom
#           (Eukaryota / Bacteria / Archaea / Viruses)
#   ERROR   Species binomial does not match the s__ rank (mislabel / cross-organism row)
#   ERROR   contrasting clades: the same taxon placed under different parents
#           across custom rows (internal inconsistency; corrupts the taxdump)
#   ERROR   lineage conflicts with the reference taxonomy: a known taxon placed
#           under any parent it does not have there (only with a reference)
#   WARNING new coherent lineage: a lineage that does not conflict but introduces
#           taxa the reference does not know - listed so novel additions can be
#           reviewed (the only warning; only with a reference)
#
# Custom proteomes are not required to be eukaryotic. Each row is checked against
# the reference matching its own domain, built WITHOUT the custom proteomes
# (mnemo<TAB>lineage): Eukaryota -> eukaryotes_taxonomy_ref.tsv,
# Bacteria/Archaea -> gtdb_taxonomy.tsv, Viruses -> virus_taxonomy.tsv. A domain
# with no reference provided is validated for schema only.
#
# Usage:
#   Rscript workflow/scripts/check_custom_proteomes.R resources/custom_genomes_repdb.csv
#   Rscript workflow/scripts/check_custom_proteomes.R <file> \
#     --reference euk_ref.tsv --reference-prok gtdb_tax.tsv \
#     --reference-virus virus_tax.tsv --check-fasta
#
# It can also be called from Snakemake (snakemake@input[["table"]],
# snakemake@input[["reference"]]).

suppressMessages(library(tidyverse))

REQUIRED_COLUMNS <- c("ID", "Species", "Data_type", "Fasta", "Lineage")
RESERVED_PREFIXES <- c("UP", "EP", "P10")
RANK_PREFIXES <- c("d__", "p__", "c__", "o__", "f__", "g__", "s__")
RANK_LABELS <- c("domain", "phylum", "class", "order", "family", "genus", "species")
RANK_KEYS <- c("d", "p", "c", "o", "f", "g", "s")
# Custom proteomes need not be eukaryotic: the user may add bacterial/archaeal
# (routed against GTDB) or viral (routed against ICTV/NCBI) proteomes. Each row
# is checked against the reference matching its own domain.
VALID_DOMAINS <- c("Eukaryota", "Bacteria", "Archaea", "Viruses")

# ---- reporting helpers -------------------------------------------------------
errors <- character(0)
warnings <- character(0)
add_error <- function(msg) errors[[length(errors) + 1]] <<- msg
add_warning <- function(msg) warnings[[length(warnings) + 1]] <<- msg

# first two tokens (genus species), lowercased and stripped of punctuation
# (so "Skoliomonas sp. GEMRC" and "Skoliomonas sp GEMRC" compare equal)
binomial <- function(name) {
  clean <- gsub("[^[:alnum:] ]", " ", tolower(trimws(name)))
  paste(head(strsplit(clean, "\\s+")[[1]], 2), collapse = " ")
}

# validate one lineage string; returns the 7 rank values or NULL if malformed
check_lineage <- function(lineage, where) {
  ranks <- str_split(lineage, ";")[[1]]
  if (length(ranks) != length(RANK_PREFIXES)) {
    add_error(sprintf(
      "%s: lineage has %d rank(s), expected %d (d__;p__;c__;o__;f__;g__;s__): '%s'",
      where, length(ranks), length(RANK_PREFIXES), lineage))
    return(NULL)
  }
  values <- character(length(RANK_PREFIXES))
  ok <- TRUE
  for (i in seq_along(RANK_PREFIXES)) {
    rank <- trimws(ranks[[i]])
    prefix <- RANK_PREFIXES[[i]]
    if (!startsWith(rank, prefix)) {
      add_error(sprintf("%s: expected rank prefix '%s' but found '%s'",
                        where, prefix, rank))
      ok <- FALSE
      next
    }
    value <- trimws(substring(rank, nchar(prefix) + 1))
    if (value == "") {
      add_error(sprintf("%s: rank '%s' has no value", where, prefix))
      ok <- FALSE
    }
    values[[i]] <- value
  }
  if (ok) values else NULL
}

# Build, from a reference eukaryotic taxonomy (mnemo<TAB>lineage, no custom
# entries): for each rank, the map child -> parent name(s), and the set of all
# taxon names seen at that rank. Used to check custom lineages against it.
parse_reference <- function(path) {
  ref <- suppressWarnings(suppressMessages(
    read_tsv(path, col_names = c("mnemo", "lineage"),
             col_types = cols(.default = "c"), name_repair = "minimal")))
  ref <- ref %>%
    separate(lineage, RANK_KEYS, sep = ";", fill = "right", extra = "drop") %>%
    mutate(across(all_of(RANK_KEYS), ~ trimws(gsub("^[a-z]__", "", .))))
  parents <- vector("list", length(RANK_KEYS))
  known <- vector("list", length(RANK_KEYS))
  for (i in seq_along(RANK_KEYS)) {
    col <- ref[[RANK_KEYS[i]]]
    known[[i]] <- unique(col[!is.na(col) & col != ""])
    if (i >= 2) {
      d <- ref[, c(RANK_KEYS[i], RANK_KEYS[i - 1])]
      names(d) <- c("child", "parent")
      d <- unique(d[!is.na(d$child) & d$child != "" &
                    !is.na(d$parent) & d$parent != "", ])
      parents[[i]] <- split(d$parent, d$child)
    }
  }
  list(parents = parents, known = known)
}

# Compare one custom lineage (7 rank values, prefixes stripped) against the
# reference. Any known taxon placed under a parent it does not have in the
# reference is a conflict -> ERROR. A lineage that has no such conflict but
# introduces taxa the reference does not know is a "new coherent lineage" -> the
# single WARNING, listing it so novel additions can be eyeballed.
check_against_reference <- function(values, where, lineage, ref) {
  conflict <- FALSE
  new_taxa <- character(0)
  for (i in 2:length(values)) {
    child <- values[[i]]
    parent <- values[[i - 1]]
    if (child == "") next
    if (!(child %in% ref$known[[i]])) {
      new_taxa <- c(new_taxa, sprintf("%s '%s'", RANK_LABELS[i], child))
      next
    }
    ref_parents <- ref$parents[[i]][[child]]
    if (parent != "" && !is.null(ref_parents) && !(parent %in% ref_parents)) {
      add_error(sprintf(
        "%s: %s '%s' is placed under %s '%s', but the reference taxonomy places it under '%s'",
        where, RANK_LABELS[i], child, RANK_LABELS[i - 1], parent,
        paste(ref_parents, collapse = "' / '")))
      conflict <- TRUE
    }
  }
  if (!conflict && length(new_taxa) > 0) {
    add_warning(sprintf("%s: new coherent lineage - adds %s | %s",
                        where, paste(new_taxa, collapse = ", "), lineage))
  }
}

# Check the custom table for *internal* contradictions: the same taxon name
# placed under different parents across custom rows (e.g. class 'Provora' under
# phylum 'Diaphoretickes' in one row and 'Metamonada' in another). This corrupts
# the taxdump and needs no reference, so it always runs.
check_internal_conflicts <- function(rows) {
  if (length(rows) < 2) return(invisible())
  for (i in 2:length(RANK_LABELS)) {
    m <- list()  # child -> named list(parent -> first "where" seen)
    for (r in rows) {
      child <- r$values[[i]]
      parent <- r$values[[i - 1]]
      if (child == "" || parent == "") next
      if (is.null(m[[child]])) m[[child]] <- list()
      if (is.null(m[[child]][[parent]])) m[[child]][[parent]] <- r$where
    }
    for (child in names(m)) {
      parents <- names(m[[child]])
      if (length(parents) > 1) {
        detail <- paste(sprintf("'%s' (%s)", parents, unlist(m[[child]])),
                        collapse = " vs ")
        add_error(sprintf(
          "custom table places %s '%s' under conflicting %s: %s",
          RANK_LABELS[i], child, RANK_LABELS[i - 1], detail))
      }
    }
  }
}

# ---- argument handling -------------------------------------------------------
# out_path is a marker written on success (only when run from Snakemake).
# reference_path is an optional reference eukaryotic taxonomy for the conflict
# check; when absent, only the schema checks run.
# reference paths per domain: euk (Eukaryota), prok (Bacteria + Archaea, GTDB),
# virus (Viruses, ICTV/NCBI). Any may be absent -> that domain is schema-only.
reference_paths <- list(euk = NULL, prok = NULL, virus = NULL)

if (exists("snakemake")) {
  table_path <- snakemake@input[["table"]]
  check_fasta <- isTRUE(snakemake@params[["check_fasta"]])
  out_path <- snakemake@output[[1]]
  opt_ref <- function(name) {
    r <- tryCatch(snakemake@input[[name]], error = function(e) NULL)
    if (is.null(r) || length(r) == 0 || !nzchar(r[[1]])) NULL else r[[1]]
  }
  reference_paths$euk <- opt_ref("reference")
  reference_paths$prok <- opt_ref("reference_prok")
  reference_paths$virus <- opt_ref("reference_virus")
} else {
  args <- commandArgs(trailingOnly = TRUE)
  check_fasta <- "--check-fasta" %in% args
  args <- args[args != "--check-fasta"]
  take_flag <- function(args, flag) {
    i <- which(args == flag)
    if (length(i) == 1 && length(args) >= i + 1) {
      list(value = args[[i + 1]], args = args[-c(i, i + 1)])
    } else {
      list(value = NULL, args = args)
    }
  }
  for (flag in c("--reference", "--reference-prok", "--reference-virus")) {
    r <- take_flag(args, flag); args <- r$args
    key <- c("--reference" = "euk", "--reference-prok" = "prok",
             "--reference-virus" = "virus")[[flag]]
    reference_paths[[key]] <- r$value
  }
  if (length(args) != 1) {
    stop(paste("usage: check_custom_proteomes.R <custom_proteome.tsv>",
               "[--reference <euk_tax.tsv>] [--reference-prok <gtdb_tax.tsv>]",
               "[--reference-virus <virus_tax.tsv>] [--check-fasta]"))
  }
  table_path <- args[[1]]
  out_path <- NULL
}

if (!file.exists(table_path)) {
  stop(sprintf("file not found: %s", table_path))
}
message(sprintf("Validating custom proteome table: %s\n", table_path))

# read everything as character; ignore any trailing empty header column
df <- suppressWarnings(suppressMessages(
  read_tsv(table_path, col_types = cols(.default = "c"), name_repair = "minimal")))
df <- df[, !is.na(names(df)) & names(df) != "", drop = FALSE]

# optional reference taxonomies for the conflict check, mapped to the domains
# they cover (prokaryote reference serves both Bacteria and Archaea).
ref_by_domain <- list()
usable_ref <- function(path) {
  !is.null(path) && file.exists(path) && file.info(path)$size > 0
}
for (spec in list(list("euk", "Eukaryota"),
                  list("prok", c("Bacteria", "Archaea")),
                  list("virus", "Viruses"))) {
  path <- reference_paths[[spec[[1]]]]
  if (usable_ref(path)) {
    message(sprintf("Checking %s lineages against reference: %s",
                    paste(spec[[2]], collapse = "/"), path))
    info <- parse_reference(path)
    for (d in spec[[2]]) ref_by_domain[[d]] <- info
  }
}

missing <- setdiff(REQUIRED_COLUMNS, names(df))
report_and_exit <- function(n_rows) {
  for (w in warnings) cat(sprintf("  WARNING: %s\n", w))
  for (e in errors) cat(sprintf("  ERROR:   %s\n", e))
  cat("\n")
  cat(sprintf("Checked %d custom proteome(s): %d error(s), %d warning(s).\n",
              n_rows, length(errors), length(warnings)))
  if (length(errors) > 0) {
    cat("FAILED - fix the errors above before running the pipeline.\n")
    quit(status = 1)
  }
  cat("OK - the custom proteome table is valid.\n")
  if (!is.null(out_path)) {
    writeLines(sprintf(
      "OK - %d custom proteome(s) valid, %d warning(s). Validated %s",
      n_rows, length(warnings), table_path), out_path)
  }
  quit(status = 0)
}

if (length(missing) > 0) {
  for (c in missing) add_error(sprintf("missing required column: '%s'", c))
  report_and_exit(0)
}

# ---- per-row checks ----------------------------------------------------------
seen_ids <- list()
seen_fasta <- list()
custom_lineages <- list()  # well-formed lineages, for the internal-conflict check
n_rows <- 0L

for (i in seq_len(nrow(df))) {
  row <- df[i, ]
  # data row i is file line i + 1 (header is line 1)
  lineno <- i + 1
  get <- function(name) {
    v <- row[[name]]
    if (is.na(v)) "" else trimws(v)
  }
  rid <- get("ID")
  if (rid == "" && get("Lineage") == "" && get("Fasta") == "") next  # blank line
  n_rows <- n_rows + 1L
  where <- sprintf("line %d (ID=%s)", lineno, if (rid == "") "<empty>" else rid)

  # ID: present, CUS-prefixed, unique
  if (rid == "") {
    add_error(sprintf("%s: empty ID", where))
  } else {
    if (!startsWith(rid, "CUS")) {
      add_error(sprintf(
        "%s: ID must start with 'CUS' (non-CUS IDs are silently dropped by filter_tax.R)",
        where))
    } else if (any(startsWith(rid, RESERVED_PREFIXES))) {
      add_error(sprintf("%s: ID collides with a reserved source prefix (%s)",
                        where, paste(RESERVED_PREFIXES, collapse = ", ")))
    }
    if (!is.null(seen_ids[[rid]])) {
      add_error(sprintf("%s: duplicate ID, first seen at line %d", where, seen_ids[[rid]]))
    } else {
      seen_ids[[rid]] <- lineno
    }
  }

  # Lineage: exactly 7 well-formed ranks
  lineage <- get("Lineage")
  values <- NULL
  if (lineage == "") {
    add_error(sprintf("%s: empty Lineage", where))
  } else {
    values <- check_lineage(lineage, where)
  }

  # checks on a well-formed lineage
  if (!is.null(values)) {
    domain <- values[[1]]
    species_rank <- values[[length(values)]]
    if (!(domain %in% VALID_DOMAINS)) {
      add_error(sprintf(
        "%s: domain '%s' is not a recognized superkingdom (%s)",
        where, domain, paste(VALID_DOMAINS, collapse = ", ")))
    }
    species <- get("Species")
    if (species != "" && binomial(species) != "" &&
        binomial(species) != binomial(species_rank)) {
      add_error(sprintf(
        "%s: Species '%s' does not match the s__ rank '%s' (possible mislabel or cross-organism row)",
        where, species, species_rank))
    }
    # compare against the reference taxonomy for this row's domain (if provided)
    ref_info <- ref_by_domain[[domain]]
    if (!is.null(ref_info)) {
      check_against_reference(values, where, lineage, ref_info)
    }
    custom_lineages[[length(custom_lineages) + 1]] <- list(where = where, values = values)
  }

  # Fasta path
  fasta <- get("Fasta")
  if (fasta == "") {
    add_error(sprintf("%s: empty Fasta path", where))
  } else {
    if (!is.null(seen_fasta[[fasta]])) {
      add_error(sprintf("%s: duplicate Fasta path, first seen at line %d",
                        where, seen_fasta[[fasta]]))
    } else {
      seen_fasta[[fasta]] <- lineno
    }
    if (check_fasta && (!file.exists(fasta) || file.info(fasta)$size == 0)) {
      add_error(sprintf("%s: Fasta path missing or empty: %s", where, fasta))
    }
  }
}

# internal contradictions across the custom rows (contrasting clades)
check_internal_conflicts(custom_lineages)

report_and_exit(n_rows)
