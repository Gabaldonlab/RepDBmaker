---
title: "RepDB vs nr: search cost, coverage and taxonomic resolution"
date: "08-Sep-2026"
author: "Giacomo Mutti"
output:
    rmdformats::html_clean:
        number_sections: false
        code_folding: hide
        self_contained: true
        fig_caption: true
        gallery: true
        toc_depth: 3
        highlight: kate
params:
    results: "results"
    resources: "resources"
    figdir: "results/plots"
    reference: "resources/merged.fasta.transdecoder_reduced.tsv"
    reference_name: "Cohen2024"
    reference_long: "Cohen et al. (2024)"
---

<!--
Render from the repository root:

  Rscript -e 'root <- getwd()
  dir.create(file.path(root, "results/qc"), recursive = TRUE, showWarnings = FALSE)
  rmarkdown::render("workflow/notebooks/comparison.Rmd",
                    knit_root_dir = root,
                    output_dir    = file.path(root, "results/qc"),
                    output_file   = "comparison.html")'

(Bind the root to a variable first: render() setwd()s into the notebook's
directory before it forces `output_file`.)
-->




``` r
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
# free() is used to lay out Fig. 5 and only exists from patchwork 1.2.0
stopifnot(packageVersion("patchwork") >= "1.2.0")
theme_set(theme_classic(base_size = 9))

.src <- function(p) source(if (file.exists(p)) p else file.path("..", "..", p))
.src("workflow/scripts/functions.R")
.src("workflow/scripts/palettes.R")
```


``` r
P <- as.list(params)
res <- function(...) rp(file.path(P$results, ...))
rsc <- function(...) rp(file.path(P$resources, ...))
have <- function(p) !is.null(p) && length(p) == 1L && nzchar(p) && file.exists(rp(p))

figdir <- rp(P$figdir)
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
fig <- function(name) file.path(figdir, name)

REF <- P$reference_name       # short, for axes and legends
REF_LONG <- P$reference_long  # full, for prose

# The reference series is not one of the four benchmarked databases, so it gets
# its own hue. The previous brown (#8C510A) collapsed against repdb green under
# protanopia (dE 5.4, below the 8 floor); this rose clears every pair at
# dE>=10.1 simulated and dE>=20.6 unsimulated. Re-check with
# `Rscript workflow/scripts/check_palettes.R` if you change it.
color_ref <- c("#CD5D8E")
names(color_ref) <- REF
db_colors <- c(color_benchdbs, color_ref)
DB_LEVELS <- names(db_colors)

# One identical colour scale on every panel. patchwork only merges guides whose
# scale AND key glyph match, so the levels are fixed with drop = FALSE and every
# geom draws a point key - otherwise panels A, B and D each emit their own
# near-identical legend.
scale_db <- function() {
  # patchwork keeps the direction a guide was built with and ignores
  # legend.direction in the assembly theme, so the orientation is set here
  scale_color_manual(values = db_colors, limits = DB_LEVELS, drop = FALSE,
                     na.translate = FALSE,
                     guide = guide_legend(direction = "vertical"))
}

# Panel C: green agrees, orange disagrees, and the two "no answer" outcomes are
# neutral greys (light = homolog but no LCA, dark = no homolog). Defined in
# palettes.R and verified by check_palettes.R - all six pairs clear dE 8 under
# protanopia and deuteranopia and dE 15 unsimulated.
OUTCOMES <- c(sprintf("Agrees with %s", REF), "Assigned, differs",
              "Homolog, no LCA", "No homolog")
outcome_colors <- setNames(outcome_colors_ref, OUTCOMES)

N <- list()

# panels carry their own theme rather than relying on theme_set(): the assembled
# object is written to an .Rds that eda_report.Rmd re-reads, and a ggplot resolves
# the global theme at draw time, so a saved plot would otherwise pick up
# whatever theme happens to be active in the reading session.
# base_size 8 at 7.5x5in put the tick labels near 6pt once the PDF was placed
# in a two-column layout, which is below what most journals accept. The figure
# is now sized to a full text width (7.2in) and given the height the 2x2 grid
# actually wants, so the same panels carry 10pt type.
thm <- theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    legend.key.size = unit(4.5, "mm"),
    legend.text = element_text(size = 9),
    legend.title = element_blank(),
    axis.text = element_text(size = 9),
    plot.tag = element_text(size = 12, face = "bold"),
    # panels A and B are facetted on free scales, so both facets run their axis
    # out to a wide tick; at 9pt the last label of the left facet and the first
    # of the right one collide unless the strips are held apart
    panel.spacing = unit(4, "mm"),
    plot.margin = margin(3, 3, 3, 3)
  )

show_tbl <- function(df, caption = NULL, digits = 1) {
  k <- knitr::kable(df, caption = caption, digits = digits,
                    format.args = list(big.mark = ","))
  if (requireNamespace("kableExtra", quietly = TRUE)) {
    k <- kableExtra::kable_styling(k, bootstrap_options = c("striped", "condensed"),
                                   full_width = FALSE)
  }
  k
}
fmt_n <- function(x) formatC(round(x), format = "d", big.mark = ",")
fmt_p <- function(x, d = 1) paste0(formatC(x, format = "f", digits = d), "%")
fmt_d <- function(x, d = 2) formatC(x, format = "f", digits = d)
fmt_x <- function(x, d = 2) paste0(formatC(x, format = "f", digits = d), "x")
```

# Input validation

Every per-transcript result must refer to the *current* query set. This gate
exists because `mmseqs easy-taxonomy` writes four files per run and it is easy to
end up with a `_lca.tsv` left over from an earlier sample while the `_report`
next to it is current — the two disagree silently and every downstream
"unclassified" number is then wrong.


``` r
fa <- readLines(res("comparison", "input.fa"), warn = FALSE)
query_ids <- sub(" .*", "", sub("^>", "", fa[startsWith(fa, ">")]))
N$n_queries <- length(query_ids)

# per-transcript LCA assignments
lca_files <- list.files(res("comparison", "tax"), pattern = "_lca\\.tsv$", full.names = TRUE)
lca_all <- map_dfr(lca_files, ~ read_delim(.x,
  col_names = c("seq", "taxid", "rank", "classification"),
  show_col_types = FALSE, progress = FALSE
) %>% mutate(db = gsub("_lca.tsv", "", basename(.x))))

# kraken-style aggregate reports (one row per taxon, depth-indented)
report_files <- list.files(res("comparison", "tax"), pattern = "_report$", full.names = TRUE)
report_files <- report_files[!str_detect(report_files, "tophit")]
reports <- map_dfr(report_files, ~ read_delim(.x, delim = "\t",
  col_names = c("perc", "n", "n_direct", "rank", "taxid", "nm"),
  show_col_types = FALSE, progress = FALSE
) %>% mutate(db = gsub("_report", "", basename(.x))))

validation <- bind_rows(
  lca_all %>% group_by(db) %>%
    summarise(File = "_lca.tsv", Rows = n(),
              `Matching current input` = sum(seq %in% query_ids), .groups = "drop"),
  reports %>% group_by(db) %>%
    summarise(File = "_report", Rows = sum(n_direct),
              `Matching current input` = NA_integer_, .groups = "drop")
) %>%
  mutate(
    `% current` = `Matching current input` / Rows * 100,
    Status = case_when(
      File == "_report" & Rows == N$n_queries ~ "OK",
      File == "_report" ~ "STALE (row total != query count)",
      `% current` > 99 ~ "OK",
      TRUE ~ "STALE - re-run mmseqs for this db"
    )
  ) %>%
  select(db, File, Rows, `Matching current input`, `% current`, Status) %>%
  arrange(File, db)

# databases whose per-transcript file can be trusted
LCA_DBS <- lca_all %>%
  group_by(db) %>%
  summarise(ok = mean(seq %in% query_ids) > 0.99) %>%
  filter(ok) %>% pull(db)
REPORT_DBS <- reports %>%
  group_by(db) %>% summarise(ok = sum(n_direct) == N$n_queries) %>%
  filter(ok) %>% pull(db)

lca <- filter(lca_all, db %in% LCA_DBS)
FULL_LCA <- setequal(LCA_DBS, unique(lca_all$db))

show_tbl(validation, sprintf("Consistency of each result file with the %s queries in `input.fa`.", fmt_n(N$n_queries)), digits = 1)
```



Table: Consistency of each result file with the 20,000 queries in `input.fa`.

|db         |File     |   Rows| Matching current input| % current|Status |
|:----------|:--------|------:|----------------------:|---------:|:------|
|clustnr    |_lca.tsv | 20,000|                 20,000|       100|OK     |
|clustrepdb |_lca.tsv | 20,000|                 20,000|       100|OK     |
|nr         |_lca.tsv | 20,000|                 20,000|       100|OK     |
|repdb      |_lca.tsv | 20,000|                 20,000|       100|OK     |
|clustnr    |_report  | 20,000|                     NA|        NA|OK     |
|clustrepdb |_report  | 20,000|                     NA|        NA|OK     |
|nr         |_report  | 20,000|                     NA|        NA|OK     |
|repdb      |_report  | 20,000|                     NA|        NA|OK     |




``` r
# DIAMOND hits: only the columns needed, the files are ~150-190 MB each
hits <- map_dfr(
  list.files(res("comparison", "hits"), pattern = "matches.tsv", full.names = TRUE),
  function(p) {
    read_delim(p, col_names = FALSE, col_select = c(1, 3, 13),
               show_col_types = FALSE, progress = FALSE) %>%
      set_names(c("qseqid", "pident", "qcovhsp")) %>%
      mutate(db = gsub("_matches.tsv", "", basename(p)))
  }
) %>% mutate(db = factor(db, levels = DB_LEVELS))

# reference annotation of the same transcripts
ref <- read_delim(rp(P$reference), show_col_types = FALSE, progress = FALSE) %>%
  transmute(
    seq = transcript_name,
    ref_rank = classification_level,
    ref_taxon = classification,
    ref_lineage = full_classification,
    ambiguous
  ) %>%
  distinct(seq, .keep_all = TRUE) %>%
  mutate(
    ref_domain     = str_split_i(ref_lineage, ";\\s*", 1),
    ref_supergroup = str_split_i(ref_lineage, ";\\s*", 2),
    ref_division   = str_split_i(ref_lineage, ";\\s*", 3)
  )

N$n_ref_total <- nrow(ref)
N$n_ref_unclassified <- sum(is.na(ref$ref_taxon) | ref$ref_taxon == "unclassified")
N$n_ref_ambiguous <- sum(ref$ambiguous == 1, na.rm = TRUE)
```

The query set is **20,000 predicted proteins** randomly sampled
from the western North Atlantic metatranscriptome of Cohen et al. (2024).

Throughout, the comparator labelled **Cohen2024** is *the taxonomic annotation
published by that study*, not a database we queried. It is often referred to as
"PhyloDB", but the annotation was produced with a PhyloDB-derived reference set
augmented with additional references and RefSeq, and its own LCA procedure and
thresholds. It is therefore a different *pipeline over a different reference*,
which is exactly why the raw unclassified counts are not comparable and why the
resolution profile below is the informative comparison.

# Fig. 5

## A - Search cost


``` r
bench <- list.files(res("benchmarks", "comparison"), full.names = TRUE, pattern = "\\.txt$") %>%
  map_dfr(~ read_delim(.x, delim = "\t", show_col_types = FALSE, progress = FALSE) %>%
    mutate(tag = str_remove(basename(.x), "\\.txt$"))) %>%
  separate(tag, c("db", "tool"), sep = "_(?=[^_]+$)") %>%
  mutate(
    tool = recode(tool, dmnd = "DIAMOND", mmseqs = "MMseqs2"),
    mem_gb = max_rss / 1024,
    minutes = s / 60,
    db_source = gsub("^clust", "", db),
    clustered = str_starts(db, "clust"),
    db = factor(db, levels = DB_LEVELS)
  )

cost <- bench %>%
  select(tool, db, minutes, mem_gb) %>%
  arrange(tool, db)

show_tbl(
  cost %>% rename(Tool = tool, DB = db, `Runtime (min)` = minutes, `Peak RSS (GB)` = mem_gb),
  "Wall time and peak memory of the homology search / taxonomic assignment on the full query set."
)
```



Table: Wall time and peak memory of the homology search / taxonomic assignment on the full query set.

|Tool    |DB         | Runtime (min)| Peak RSS (GB)|
|:-------|:----------|-------------:|-------------:|
|DIAMOND |nr         |          42.8|          14.0|
|DIAMOND |clustnr    |          16.1|          10.4|
|DIAMOND |repdb      |          22.7|           8.7|
|DIAMOND |clustrepdb |          18.3|           7.7|
|MMseqs2 |nr         |         103.6|         107.2|
|MMseqs2 |clustnr    |          34.1|         108.7|
|MMseqs2 |repdb      |          45.6|         108.7|
|MMseqs2 |clustrepdb |          25.4|          77.3|

``` r
ratio <- cost %>%
  filter(db %in% c("nr", "repdb")) %>%
  pivot_wider(names_from = db, values_from = c(minutes, mem_gb)) %>%
  mutate(`Runtime speedup` = minutes_nr / minutes_repdb,
         `Memory reduction` = mem_gb_nr / mem_gb_repdb)
ratio_c <- cost %>%
  filter(db %in% c("clustnr", "clustrepdb")) %>%
  pivot_wider(names_from = db, values_from = c(minutes, mem_gb)) %>%
  mutate(`Runtime speedup` = minutes_clustnr / minutes_clustrepdb,
         `Memory reduction` = mem_gb_clustnr / mem_gb_clustrepdb)

gv <- function(df, tool, col) df[[col]][df$tool == tool]
N$speed_dmnd  <- gv(ratio, "DIAMOND", "Runtime speedup")
N$speed_mmseqs <- gv(ratio, "MMseqs2", "Runtime speedup")
N$mem_dmnd    <- gv(ratio, "DIAMOND", "Memory reduction")
N$mem_mmseqs  <- gv(ratio, "MMseqs2", "Memory reduction")
N$mem_mmseqs_nr    <- gv(filter(cost, db == "nr"), "MMseqs2", "mem_gb")
N$mem_mmseqs_repdb <- gv(filter(cost, db == "repdb"), "MMseqs2", "mem_gb")

show_tbl(
  bind_rows(
    ratio %>% transmute(Comparison = "RepDB vs nr", Tool = tool, `Runtime speedup`, `Memory reduction`),
    ratio_c %>% transmute(Comparison = "ClusteredRepDB vs ClusteredNR", Tool = tool, `Runtime speedup`, `Memory reduction`)
  ),
  "Ratios (>1 means RepDB is faster / leaner).", digits = 2
)
```



Table: Ratios (>1 means RepDB is faster / leaner).

|Comparison                    |Tool    | Runtime speedup| Memory reduction|
|:-----------------------------|:-------|---------------:|----------------:|
|RepDB vs nr                   |DIAMOND |            1.89|             1.61|
|RepDB vs nr                   |MMseqs2 |            2.27|             0.99|
|ClusteredRepDB vs ClusteredNR |DIAMOND |            0.88|             1.35|
|ClusteredRepDB vs ClusteredNR |MMseqs2 |            1.34|             1.41|

``` r
# A and B suppress their keys and D carries the single collected legend, so the
# fill scale here never has to merge with D's colour scale
plot_bench <- bench %>%
  ggplot(aes(mem_gb, minutes, fill = db, group = db_source)) +
  geom_line(color = "grey75", linewidth = .3) +
  geom_point(size = 2.5, pch=21) +
  facet_wrap(~tool, scales = "free") +
  scale_fill_manual(values = db_colors) +
  # the memory axes span very different ranges; give the outermost tick room so
  # it is not clipped by the panel edge in the narrow left column
  scale_x_continuous(n.breaks = 4, expand = expansion(mult = 0.12)) +
  labs(x = "Peak memory (GB)", y = "Runtime (min)") +
  thm + theme(legend.position = "none")   # the shared key lives on panel D
```

Against nr, RepDB is **1.89x faster for DIAMOND** and
**2.27x faster for MMseqs2**. DIAMOND memory drops
1.61x, but MMseqs2 peak memory is essentially unchanged
(109 GB vs 107 GB): the
`easy-taxonomy` prefilter sizes its working set to the available RAM rather than
to the database, so database size does not translate into a memory saving there.
Clustering is what reduces it (see the table above).

## B - Homology search results


``` r
best <- hits %>%
  group_by(db, qseqid) %>%
  slice_max(pident, n = 1, with_ties = FALSE) %>%
  ungroup()

coverage <- hits %>%
  group_by(db) %>%
  summarise(
    Hits = n(),
    `Queries with a hit` = n_distinct(qseqid),
    `% of queries` = n_distinct(qseqid) / N$n_queries * 100,
    `Hits per query` = n() / n_distinct(qseqid),
    `Median identity` = median(pident),
    `Median coverage` = median(qcovhsp),
    .groups = "drop"
  )

N$hit_repdb <- coverage$`% of queries`[coverage$db == "repdb"]
N$hit_nr    <- coverage$`% of queries`[coverage$db == "nr"]
N$nohit_repdb <- N$n_queries - coverage$`Queries with a hit`[coverage$db == "repdb"]

show_tbl(coverage, "DIAMOND homology search: how many of the queries each database can place at all.")
```



Table: DIAMOND homology search: how many of the queries each database can place at all.

|db         |      Hits| Queries with a hit| % of queries| Hits per query| Median identity| Median coverage|
|:----------|---------:|------------------:|------------:|--------------:|---------------:|---------------:|
|nr         | 1,449,985|             17,811|         89.1|           81.4|            60.5|            89.5|
|clustnr    | 1,413,589|             17,841|         89.2|           79.2|            57.4|            88.6|
|repdb      | 1,579,187|             19,243|         96.2|           82.1|            60.4|            89.2|
|clustrepdb | 1,554,358|             19,247|         96.2|           80.8|            58.4|            88.7|

``` r
plot_dmnd <- best %>%
  select(db, pident, qcovhsp) %>%
  pivot_longer(c(pident, qcovhsp)) %>%
  mutate(name = factor(recode(name, pident = "% identity", qcovhsp = "% coverage"),
                       levels = c("% identity", "% coverage"))) %>%
  ggplot(aes(value, color = db)) +
  geom_density(linewidth = .5, key_glyph = "point") +
  facet_wrap(~name, scales = "free") +
  # both facets run to 100; without room at the edge that last tick label is
  # clipped by the panel border
  scale_x_continuous(expand = expansion(mult = 0.06)) +
  scale_db() +
  labs(x = "", y = "Density") +
  thm + theme(legend.position = "none")   # the shared key lives on panel D
```

**RepDB finds a homolog for 96.2% of the queries against nr's
89.1%**, despite being roughly a third of the size. Only
757 queries have no RepDB homolog at all.


``` r
id_bands <- function(d) {
  d %>%
    mutate(band = cut(pident, c(0, 30, 50, 70, 90, 100),
                      labels = c("<30", "30-50", "50-70", "70-90", ">90"))) %>%
    count(db, band) %>%
    group_by(db) %>% mutate(pct = n / sum(n) * 100) %>% ungroup() %>%
    select(-n) %>% pivot_wider(names_from = db, values_from = pct) %>%
    select(band, any_of(DB_LEVELS))
}

# best hit per query - what Fig. 5B shows, and what matters for annotation
bands <- id_bands(best)
# all hits - where the effect of clustering on redundancy is visible
bands_all <- id_bands(hits)

pick_band <- function(tbl, db, b) sum(tbl[[db]][tbl$band %in% b])
N$lowid_repdb <- pick_band(bands, "repdb", c("<30", "30-50"))
N$lowid_nr    <- pick_band(bands, "nr",    c("<30", "30-50"))
N$midid_repdb <- pick_band(bands, "repdb", "70-90")
N$midid_nr    <- pick_band(bands, "nr",    "70-90")
for (d in DB_LEVELS[DB_LEVELS %in% names(bands_all)]) {
  N[[paste0("hi90_", d)]] <- pick_band(bands_all, d, ">90")
}

show_tbl(bands, "Identity of the **best hit per query** (% of queries with a hit). This is what Fig. 5B shows.")
```



Table: Identity of the **best hit per query** (% of queries with a hit). This is what Fig. 5B shows.

|band  |   nr| clustnr| repdb| clustrepdb|
|:-----|----:|-------:|-----:|----------:|
|<30   |  0.5|     0.5|   0.3|        0.2|
|30-50 | 18.2|    17.9|  11.4|       11.1|
|50-70 | 27.6|    27.8|  28.4|       28.5|
|70-90 | 22.1|    23.5|  28.4|       29.6|
|>90   | 31.6|    30.3|  31.6|       30.6|

``` r
show_tbl(bands_all, "Identity of **all hits** (% of hits). Clustering is visible here rather than in the best-hit view: a query's best hit is usually retained as a cluster representative, so removing near-identical redundancy barely moves the best-hit distribution but strips most of the >90% identity hits from the full output.")
```



Table: Identity of **all hits** (% of hits). Clustering is visible here rather than in the best-hit view: a query's best hit is usually retained as a cluster representative, so removing near-identical redundancy barely moves the best-hit distribution but strips most of the >90% identity hits from the full output.

|band  |   nr| clustnr| repdb| clustrepdb|
|:-----|----:|-------:|-----:|----------:|
|<30   |  2.6|     2.8|   2.0|        2.1|
|30-50 | 30.9|    33.6|  30.0|       31.8|
|50-70 | 29.9|    35.4|  33.2|       36.1|
|70-90 | 25.3|    24.5|  24.2|       24.4|
|>90   | 11.3|     3.7|  10.6|        5.5|

RepDB's best hits are *closer*: 11.6% fall below 50% identity
versus 18.8% for nr, and 28.4% lie between 70%
and 90% versus 22.1%. Query coverage, by contrast, is
indistinguishable between the two.

Across **all** hits, clustering strips the high-identity redundancy as intended:
hits above 90% identity fall from 11.3% to 3.7%
of the total in ClusteredNR, and from 10.6% to
5.5% in ClusteredRepDB.

## C - What happened to every query

Panels B and D describe the queries a database *does* place. This one accounts
for all of them, so the coverage rate, the agreement rate and the reason a
transcript ends up unlabelled can be read off a single decomposition. It is
restricted to the queries Cohen2024 itself placed in a domain, so every category
shares one denominator.


``` r
# lineages rebuilt from the depth-indented report, so an LCA taxid resolves to a
# domain without needing the taxdump
lineage_of <- function(rep) {
  rep <- rep %>% mutate(depth = (nchar(nm) - nchar(str_trim(nm, "left"))) %/% 2,
                        nm = str_trim(nm))
  stack <- character(0); out <- character(nrow(rep))
  for (i in seq_len(nrow(rep))) {
    d <- rep$depth[i]
    if (length(stack) > d) stack <- stack[seq_len(d)]
    stack <- c(stack, rep$nm[i])
    out[i] <- paste(stack, collapse = ";")
  }
  tibble(taxid = rep$taxid, lineage = out)
}
lin <- reports %>% filter(db %in% REPORT_DBS) %>%
  group_split(db) %>%
  map_dfr(~ lineage_of(.x) %>% mutate(db = .x$db[1]))

DOMAINS <- c("Bacteria", "Archaea", "Eukaryota", "Viruses")

agreement <- lca %>%
  left_join(distinct(lin, db, taxid, lineage), by = c("db", "taxid")) %>%
  mutate(db_domain = map_chr(str_split(lineage, ";"),
    ~ { h <- intersect(.x, DOMAINS); if (length(h)) h[1] else NA_character_ })) %>%
  inner_join(ref, by = "seq") %>%
  filter(ref_domain %in% DOMAINS)

outcome <- agreement %>%
  left_join(hits %>% distinct(db, qseqid) %>%
              mutate(db = as.character(db), has_hit = TRUE),
            by = c("db", "seq" = "qseqid")) %>%
  mutate(
    has_hit = replace_na(has_hit, FALSE),
    outcome = case_when(
      !is.na(db_domain) & db_domain == ref_domain ~ OUTCOMES[1],
      !is.na(db_domain)                           ~ OUTCOMES[2],
      has_hit                                     ~ OUTCOMES[3],
      TRUE                                        ~ OUTCOMES[4]
    )
  )

outcome_tbl <- outcome %>%
  count(db, outcome) %>%
  group_by(db) %>% mutate(pct = n / sum(n) * 100) %>% ungroup()

show_tbl(
  outcome_tbl %>% select(-n) %>%
    pivot_wider(names_from = outcome, values_from = pct) %>%
    rename(DB = db) %>% select(DB, any_of(OUTCOMES)),
  sprintf("What happened to each of the %s queries %s placed in a domain (%% of that set).",
          fmt_n(n_distinct(outcome$seq)), REF)
)
```



Table: What happened to each of the 19,503 queries Cohen2024 placed in a domain (% of that set).

|DB         | Agrees with Cohen2024| Assigned, differs| Homolog, no LCA| No homolog|
|:----------|---------------------:|-----------------:|---------------:|----------:|
|clustnr    |                  74.9|               3.0|            12.8|        9.3|
|clustrepdb |                  84.8|               1.9|             9.8|        3.4|
|nr         |                  74.9|               3.1|            12.6|        9.4|
|repdb      |                  84.5|               1.7|            10.3|        3.5|

``` r
plot_outcome <- outcome_tbl %>%
  mutate(db = factor(db, levels = DB_LEVELS),
         outcome = factor(outcome, levels = rev(OUTCOMES))) %>%
  ggplot(aes(db, pct, fill = outcome)) +
  geom_col(position = "stack", width = .75) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  scale_fill_manual(values = outcome_colors, breaks = OUTCOMES,
                    guide = guide_legend(direction = "vertical")) +
  labs(x = "", y = "% of queries") +
  thm +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
```

## Kingdom composition (supplementary)


``` r
ref_kingdom <- ref %>%
  mutate(nm = replace_na(ref_domain, "unclassified")) %>%
  count(nm) %>% mutate(db = REF)

plot_kings <- reports %>%
  filter(db %in% REPORT_DBS) %>%
  mutate(depth = (nchar(nm) - nchar(str_trim(nm, "left"))) %/% 2, nm = str_trim(nm)) %>%
  filter(nm %in% c("Bacteria", "Archaea", "Eukaryota", "Viruses", "unclassified")) %>%
  group_by(db, nm) %>% summarise(n = sum(n), .groups = "drop") %>%
  bind_rows(ref_kingdom) %>%
  mutate(
    db = factor(db, levels = DB_LEVELS),
    nm = factor(nm, levels = c("unclassified", names(color_kingdoms)))
  ) %>%
  ggplot(aes(db, n, fill = nm)) +
  geom_col(position = "fill", width = .75) +
  labs(x = "", y = "Proportion of queries") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  scale_fill_manual(values = c("unclassified" = "grey65", color_kingdoms)) +
  thm +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
```

## D - Taxonomic resolution

Resolution is what makes the comparison with the reference annotation
interpretable: the raw count of labelled transcripts says nothing about how
deeply they were labelled.


``` r
# the four databases use different rank vocabularies (nr carries the full NCBI
# set, RepDB a compact 7-rank taxdump, the reference a PhyloDB-derived scheme).
# Fold
# every label onto the seven canonical ranks before comparing.
CANON <- c(
  superkingdom = "domain", domain = "domain", kingdom = "kingdom",
  phylum = "phylum", subphylum = "phylum",
  superclass = "class", class = "class", subclass = "class", infraclass = "class",
  superorder = "order", order = "order", suborder = "order", infraorder = "order",
  superfamily = "family", family = "family", subfamily = "family",
  genus = "genus", subgenus = "genus",
  `species group` = "species", species = "species", subspecies = "species",
  strain = "species", forma = "species"
)
RANK_LEVELS <- c("species", "genus", "family", "order", "class", "phylum",
                 "kingdom", "domain", "unresolved")

canonise <- function(x) ifelse(x %in% names(CANON), unname(CANON[x]), "unresolved")

# n_direct in the report = transcripts assigned exactly at that taxon, so the
# aggregate report reproduces the per-transcript rank profile without needing
# the (sometimes stale) _lca.tsv
db_profile <- reports %>%
  filter(db %in% REPORT_DBS) %>%
  mutate(cr = canonise(rank)) %>%
  group_by(db, cr) %>% summarise(n = sum(n_direct), .groups = "drop")

# the reference uses a PhyloDB-derived scheme: supergroup and division have no exact
# NCBI equivalent, so they are folded onto kingdom and phylum. family / genus /
# species mean the same thing in both, which is where the comparison is safe.
ref_profile <- ref %>%
  mutate(cr = case_when(
    ref_rank %in% c("species", "genus", "family", "order", "class") ~ ref_rank,
    ref_rank == "division"   ~ "phylum",
    ref_rank == "supergroup" ~ "kingdom",
    TRUE ~ "unresolved"
  )) %>%
  count(cr) %>% mutate(db = REF)

profile <- bind_rows(db_profile, ref_profile) %>%
  mutate(cr = factor(cr, levels = RANK_LEVELS)) %>%
  group_by(db) %>% arrange(cr, .by_group = TRUE) %>%
  mutate(pct = n / sum(n) * 100, cum = cumsum(pct)) %>%
  ungroup()

res_tbl <- profile %>%
  filter(cr != "unresolved") %>%
  select(db, cr, cum) %>%
  pivot_wider(names_from = db, values_from = cum) %>%
  rename(`Rank or finer` = cr)

show_tbl(res_tbl, "Cumulative percentage of the query set assigned at each rank *or finer*. Read down: the reference labels almost everything, but very little of it deeply.")
```



Table: Cumulative percentage of the query set assigned at each rank *or finer*. Read down: the reference labels almost everything, but very little of it deeply.

|Rank or finer | Cohen2024| clustnr| clustrepdb|   nr| repdb|
|:-------------|---------:|-------:|----------:|----:|-----:|
|species       |       2.4|    39.2|       42.8| 31.4|  35.2|
|genus         |      10.3|    42.3|       51.9| 35.5|  50.7|
|family        |      26.6|    46.5|       66.9| 39.8|  66.3|
|order         |      46.7|    51.0|       70.3| 45.0|  69.7|
|class         |      72.9|    57.6|       73.3| 52.8|  73.2|
|phylum        |      81.5|    62.0|       79.5| 59.7|  79.2|
|kingdom       |      97.5|    63.5|         NA| 61.5|    NA|
|domain        |        NA|    72.5|       86.2| 72.4|  85.7|

``` r
pv <- function(d, r) {
  v <- profile$cum[profile$db == d & profile$cr == r]
  if (length(v) == 0) NA else v
}
N$fam_repdb <- pv("repdb", "family")
N$fam_nr    <- pv("nr", "family")
N$fam_ref   <- pv(REF, "family")
N$sp_repdb  <- pv("repdb", "species")
N$sp_nr     <- pv("nr", "species")
N$sp_ref    <- pv(REF, "species")
N$gen_repdb <- pv("repdb", "genus")
N$gen_ref   <- pv(REF, "genus")
N$any_repdb <- pv("repdb", "domain")
N$any_nr    <- pv("nr", "domain")
N$any_ref   <- pv(REF, "kingdom")

# Plot only species..phylum. Above phylum the schemes stop corresponding -
# the reference's "supergroup" is not an NCBI kingdom and the mmseqs databases emit
# superkingdom but never kingdom - so a line through those ranks would compare
# vocabulary, not performance. The "any label at all" figures are in the table.
PLOT_RANKS <- rev(c("species", "genus", "family", "order", "class", "phylum"))

plot_res <- profile %>%
  filter(cr %in% PLOT_RANKS) %>%
  mutate(
    cr = factor(cr, levels = PLOT_RANKS),
    db = factor(db, levels = DB_LEVELS)
  ) %>%
  ggplot(aes(as.integer(cr), cum, color = db, group = db)) +
  geom_line(linewidth = .6, key_glyph = "point") +
  geom_point(size = 1.5) +
  scale_x_continuous(
    breaks = seq_along(PLOT_RANKS), labels = str_to_title(PLOT_RANKS),
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  scale_db() +
  labs(x = "Assigned at rank or finer", y = "% of queries") +
  thm +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
```

At **family level or finer** — the deepest rank where the reference's and the
NCBI/UniEuk vocabularies mean the same thing — RepDB places **66.3%** of the
queries, nr 39.8%, and the reference annotation only
26.6%.

## Assembled figure


``` r
# same shape as the original figure - the two cost/quality panels stacked on the
# left, the taxonomy panels to the right - with A and B given less width now that
# there are four panels. A and B carry no key of their own; the collected legends
# (kingdoms from C, databases from D) sit at the right of the whole assembly.
#
# Built as ONE flat 2x2 grid, not as `(A / B) | (C / D)`: nesting makes two
# independent column groups, patchwork sizes their rows separately, and the D
# tag ends up lower than B. wrap_plots fills row-major, so the plots go in
# reading order and tag_levels spells out the letters they should carry.
#
# free(..., "space", "b") lets panel A keep its own bottom furniture. In a flat
# grid the axis-label row is shared across a row, so C's angled tick labels
# would otherwise size it and strand A's "Peak memory (GB)" title an inch below
# its own ticks.
#
# Geometry: 7.5in wide, which is what panel B needs for the outermost "100" of
# both facets to clear the panel border. The height is solved, not guessed - the
# panels are measured off the rendered PNG and 6.19in is what reproduces the
# panel height of the earlier 7.2x5.94 draft, the flat grid having given ~6% of
# it back to the now-shared axis rows. The original 7.5x5 was close to 3:2, which
# squeezed each panel to ~2.4in tall and forced 8pt type; the panels are now
# nearly square, which is also the right shape for D (a rank profile read left to
# right against a full 0-100 axis).
plot_comp <- wrap_plots(
    free(plot_bench, type = "space", side = "b"), plot_outcome,
    plot_dmnd, plot_res,
    ncol = 2
  ) +
  plot_layout(widths = c(1.3, 1), guides = "collect") +
  # Both collected keys sit in a column on the right. legend.spacing.y is the
  # only lever ggplot gives for spreading stacked keys - there is no
  # space-between - and 40mm is tuned so the outcome key lands beside panel C
  # and the database key beside panel D rather than both bunching at the
  # vertical centre. It is tied to the 6.19in height below; re-tune it if that
  # changes.
  #
  # This has to go through plot_annotation(theme = ), NOT `& theme(...)`: `&`
  # applies to every panel and would overwrite the legend.position = "none" that
  # A and B set, so the figure would come back with three redundant database
  # keys instead of one.
  plot_annotation(
    # wrap_plots filled row-major (bench, outcome, dmnd, res), so the tags are
    # given explicitly rather than left to run A-B-C-D down that order
    tag_levels = list(c("A", "C", "B", "D")),
    theme = theme(legend.position = "right", legend.box = "vertical",
                  legend.box.just = "left", legend.margin = margin(0, 0, 0, 0),
                  legend.box.spacing = unit(3, "mm"),
                  legend.spacing.y = unit(45, "mm"))
  )

saveRDS(plot_comp, fig("comparison.Rds"))
ggsave(fig("Fig5.pdf"), plot_comp, width = 7.5, height = 6.19)
plot_comp
```

<img src="/Users/gmutti/Desktop/projects/RepDBmaker/results/qc/comparison_files/figure-html/fig5-1.png" alt="" width="720" />

*A) Wall time and peak memory of the search on 20,000 predicted
proteins; the grey line joins each database to its clustered version. B) Identity
and query coverage of the best DIAMOND hit per query. C) Outcome of every query
Cohen2024 placed in a domain: whether the database agreed with the reference
at domain level, assigned a different domain, found a homolog without resolving
an LCA, or found no homolog at all. D) Cumulative percentage of queries assigned at each rank or finer,
over the range where the reference, NCBI and UniEuk/GTDB rank vocabularies
correspond.*

# Why the reference leaves fewer transcripts unclassified

The reference annotation labels more transcripts than RepDB does, which invites
the reading that RepDB's downsampling has removed sequences needed to annotate
environmental data. It has not: RepDB is ~3.5x smaller than nr and still finds
a homolog for more of the queries (96.2% against
89.1%, table in section B above). The difference is in what
"classified" means on each side, and the sections below take that apart.

## 1. What "unclassified" actually means

`unclassified` conflates two very different failures: *no homolog exists in the
database* (a genuine coverage gap) and *homologs exist but the LCA declined to
commit to a taxon* (a property of the LCA over a taxonomically broad database).
Splitting them separates a genuine coverage gap from an algorithmic artifact.


``` r
UNINFORMATIVE <- c("unclassified", "root", "cellular organisms")

has_hit <- hits %>% distinct(db, qseqid) %>%
  mutate(db = as.character(db), has_hit = TRUE)

reason <- lca %>%
  mutate(unresolved = is.na(classification) | classification %in% UNINFORMATIVE) %>%
  left_join(has_hit, by = c("db", "seq" = "qseqid")) %>%
  mutate(has_hit = replace_na(has_hit, FALSE))

reason_tbl <- reason %>%
  filter(unresolved) %>%
  count(db, has_hit) %>%
  group_by(db) %>%
  mutate(`% of unresolved` = n / sum(n) * 100) %>%
  ungroup() %>%
  transmute(db,
    Reason = if_else(has_hit, "Homolog found, LCA did not commit", "No homolog in the database"),
    Transcripts = n, `% of unresolved`
  )

N$lca_artifact_pct <- reason_tbl$`% of unresolved`[
  reason_tbl$db == "repdb" & str_starts(reason_tbl$Reason, "Homolog")]

show_tbl(reason_tbl, "Why a transcript ends up without a taxonomic label. DIAMOND and MMseqs2 use different thresholds, so 'homolog found' is a proxy.")
```



Table: Why a transcript ends up without a taxonomic label. DIAMOND and MMseqs2 use different thresholds, so 'homolog found' is a proxy.

|db         |Reason                            | Transcripts| % of unresolved|
|:----------|:---------------------------------|-----------:|---------------:|
|clustnr    |No homolog in the database        |       1,829|            40.6|
|clustnr    |Homolog found, LCA did not commit |       2,674|            59.4|
|clustrepdb |No homolog in the database        |         680|            24.7|
|clustrepdb |Homolog found, LCA did not commit |       2,072|            75.3|
|nr         |No homolog in the database        |       1,852|            41.4|
|nr         |Homolog found, LCA did not commit |       2,626|            58.6|
|repdb      |No homolog in the database        |         682|            23.9|
|repdb      |Homolog found, LCA did not commit |       2,174|            76.1|

For RepDB, **76.1% of unresolved transcripts do have a
homolog** — the LCA simply spans too many taxa to commit. That is a property of
the LCA over a tree-of-life database, not a missing-sequence problem.

## 2. Resolution, not coverage

The reference labels more transcripts, but far more shallowly. It comes from a
protist-focused reference set and a different LCA procedure, and its bar for
calling something "classified" is much lower than ours.


``` r
show_tbl(
  res_tbl %>% filter(`Rank or finer` %in% c("species", "genus", "family", "class")),
  "Cumulative % of queries assigned at rank or finer."
)
```



Table: Cumulative % of queries assigned at rank or finer.

|Rank or finer | Cohen2024| clustnr| clustrepdb|   nr| repdb|
|:-------------|---------:|-------:|----------:|----:|-----:|
|species       |       2.4|    39.2|       42.8| 31.4|  35.2|
|genus         |      10.3|    42.3|       51.9| 35.5|  50.7|
|family        |      26.6|    46.5|       66.9| 39.8|  66.3|
|class         |      72.9|    57.6|       73.3| 52.8|  73.2|

The reference assigns *some* label to 97.5% of transcripts against RepDB's 85.7%. But it reaches species for only 2.4% (RepDB 35.2%), genus for 10.3% (RepDB 50.7%) and family for 26.6% (RepDB 66.3%). Below family the ordering reverses and stays reversed: the reference's advantage is entirely in shallow labels.

Two further caveats on the reference labels themselves:


``` r
ref_counts <- c(N$n_ref_total, N$n_ref_ambiguous, N$n_ref_unclassified)

show_tbl(
  tibble(
    Property = c("Transcripts annotated", "Flagged `ambiguous` by the original study",
                 "Left unclassified by the original study"),
    Transcripts = ref_counts,
    `%` = ref_counts / N$n_ref_total * 100
  ),
  "Quality of the reference annotation used as ground truth."
)
```



Table: Quality of the reference annotation used as ground truth.

|Property                                  | Transcripts|     %|
|:-----------------------------------------|-----------:|-----:|
|Transcripts annotated                     |      20,000| 100.0|
|Flagged `ambiguous` by the original study |       9,208|  46.0|
|Left unclassified by the original study   |         497|   2.5|

## 3. Where RepDB and the reference disagree


``` r
agr_tbl <- agreement %>%
  group_by(db) %>%
  summarise(
    `Reference-placed queries` = n(),
    `Also placed by the DB` = sum(!is.na(db_domain)),
    Agree = sum(!is.na(db_domain) & db_domain == ref_domain),
    `% agreement where both assign` = Agree / `Also placed by the DB` * 100,
    .groups = "drop"
  )

N$domain_agreement <- agr_tbl$`% agreement where both assign`[agr_tbl$db == "repdb"]

show_tbl(agr_tbl, "Domain-level agreement with Cohen2024 (Fig. 5C). Domain is the deepest rank at which the reference, NCBI and GTDB nomenclature are directly comparable - below it, RepDB's GTDB prokaryotic names are simply different strings from the reference's NCBI names, so a string comparison would measure nomenclature, not accuracy.")
```



Table: Domain-level agreement with Cohen2024 (Fig. 5C). Domain is the deepest rank at which the reference, NCBI and GTDB nomenclature are directly comparable - below it, RepDB's GTDB prokaryotic names are simply different strings from the reference's NCBI names, so a string comparison would measure nomenclature, not accuracy.

|db         | Reference-placed queries| Also placed by the DB|  Agree| % agreement where both assign|
|:----------|------------------------:|---------------------:|------:|-----------------------------:|
|clustnr    |                   19,503|                15,183| 14,602|                          96.2|
|clustrepdb |                   19,503|                16,914| 16,534|                          97.8|
|nr         |                   19,503|                15,196| 14,600|                          96.1|
|repdb      |                   19,503|                16,813| 16,480|                          98.0|

Where both commit to a domain they agree
98.0% of the time, so
the extra transcripts RepDB resolves are not being bought with wrong calls.

## 4. Eukaryotic composition below domain level

Domain is a blunt comparison, and for a marine metatranscriptome the interesting
question is whether RepDB recovers the same *community*. Going deeper is possible,
but not by comparing taxon strings directly: UniEuk and the reference disagree on
names far more often than they disagree on biology.


``` r
euk_lineage_names <- function(x) {
  unique(unlist(str_split(x[str_detect(replace_na(x, ""), "Eukaryota")], ";\\s*")))
}
ref_vocab <- euk_lineage_names(ref$ref_lineage)

vocab_demo <- tibble(
  `Reference name` = c("Bacillariophyta", "Dinophyta", "Karlodinium micrum",
                       "Hacrobia", "Archaeplastida", "Prymnesiales",
                       "Dinophyceae_X", "Prymnesiophyceae"),
  Meaning = c("diatoms", "dinoflagellates", "a dinoflagellate species",
              "Haptophyta + Cryptophyta", "plants + red/green algae",
              "a haptophyte order", "placeholder: unassigned within Dinophyceae",
              "haptophytes")
)
```


``` r
db_vocab <- reports %>% filter(db == "repdb") %>% pull(nm) %>% str_trim() %>% unique()

show_tbl(
  vocab_demo %>%
    mutate(`Present in RepDB (UniEuk)?` = if_else(`Reference name` %in% db_vocab, "yes", "NO"),
           `UniEuk equivalent` = c("Diatomeae", "Dinoflagellata", "Karlodinium veneficum",
                                   "(not used; split)", "(not used; split)", "(not used)",
                                   "(not a taxon)", "Prymnesiophyceae")),
  "Why a raw string comparison fails below domain: synonyms, defunct groupings, junior synonyms and PR2-style placeholders. Half of these name the same biology under a different label."
)
```



Table: Why a raw string comparison fails below domain: synonyms, defunct groupings, junior synonyms and PR2-style placeholders. Half of these name the same biology under a different label.

|Reference name     |Meaning                                    |Present in RepDB (UniEuk)? |UniEuk equivalent     |
|:------------------|:------------------------------------------|:--------------------------|:---------------------|
|Bacillariophyta    |diatoms                                    |NO                         |Diatomeae             |
|Dinophyta          |dinoflagellates                            |NO                         |Dinoflagellata        |
|Karlodinium micrum |a dinoflagellate species                   |NO                         |Karlodinium veneficum |
|Hacrobia           |Haptophyta + Cryptophyta                   |NO                         |(not used; split)     |
|Archaeplastida     |plants + red/green algae                   |NO                         |(not used; split)     |
|Prymnesiales       |a haptophyte order                         |NO                         |(not used)            |
|Dinophyceae_X      |placeholder: unassigned within Dinophyceae |NO                         |(not a taxon)         |
|Prymnesiophyceae   |haptophytes                                |yes                        |Prymnesiophyceae      |

The consequence is that a per-transcript string concordance below domain would
mostly measure nomenclature. What *is* meaningful is to map both vocabularies onto
an explicit panel of marine eukaryotic groups, with the synonyms declared, and
compare composition.


``` r
# Each group lists the names that denote it in EITHER vocabulary. Within a
# lineage the first matching group wins, so entries run fine -> coarse; the
# "other X" entries therefore catch assignments that stop at the parent.
EUK_GROUPS <- list(
  "Dinophyceae"          = c("Dinophyceae", "Dinophyta", "Dinoflagellata"),
  "Ciliophora"           = c("Ciliophora"),
  "Apicomplexa+rel."     = c("Apicomplexa", "Perkinsea", "Chromerida"),
  "Haptophyta"           = c("Haptophyta", "Prymnesiophyceae", "Pavlovophyceae"),
  "Cryptophyta"          = c("Cryptophyta", "Cryptophyceae", "Cryptomonadales"),
  "Bacillariophyta"      = c("Bacillariophyta", "Bacillariophyceae", "Mediophyceae",
                             "Coscinodiscophyceae", "Diatomeae"),
  "Pelagophyceae"        = c("Pelagophyceae"),
  "Dictyochophyceae"     = c("Dictyochophyceae"),
  "other Stramenopiles"  = c("Stramenopiles", "Ochrophyta", "Oomycota", "Bicosoecida",
                             "Labyrinthulomycetes"),
  "Chlorophyta"          = c("Chlorophyta", "Mamiellophyceae", "Chlorophyceae",
                             "Trebouxiophyceae", "Chlorodendrophyceae", "Prasinodermophyta"),
  "other Archaeplastida" = c("Streptophyta", "Rhodophyta", "Glaucocystophyta",
                             "Glaucophyta", "Archaeplastida"),
  "Rhizaria"             = c("Rhizaria", "Cercozoa", "Foraminifera", "Radiolaria", "Retaria"),
  "Opisthokonta"         = c("Opisthokonta", "Metazoa", "Fungi", "Choanoflagellida",
                             "Choanoflagellatea", "Choanoflagellata", "Ascomycota",
                             "Basidiomycota", "Chytridiomycota", "Arthropoda", "Cnidaria"),
  "Amoebozoa"            = c("Amoebozoa", "Lobosa", "Conosa", "Arcellinida"),
  "Discoba/Excavata"     = c("Discoba", "Euglenozoa", "Excavata", "Heterolobosea", "Metamonada"),
  "other Alveolata"      = c("Alveolata")
)
# nesting among the panel entries, so a coarse call is not scored as a conflict
EUK_PARENT <- c(
  "Dinophyceae" = "other Alveolata", "Ciliophora" = "other Alveolata",
  "Apicomplexa+rel." = "other Alveolata",
  "Bacillariophyta" = "other Stramenopiles", "Pelagophyceae" = "other Stramenopiles",
  "Dictyochophyceae" = "other Stramenopiles",
  "Chlorophyta" = "other Archaeplastida"
)
UNRES <- "Eukaryota, unresolved"

assign_euk_group <- function(lineages) {
  parts <- str_split(lineages, ";\\s*")
  map_chr(parts, function(p) {
    for (g in names(EUK_GROUPS)) if (any(EUK_GROUPS[[g]] %in% p)) return(g)
    UNRES
  })
}
```


``` r
# A synonym that exists in neither vocabulary is dead weight; one that exists
# only in the reference makes the group look absent from RepDB. Both are silent
# failures, so the mapping is audited rather than trusted.
vocab_of <- function(d) reports %>% filter(db == d) %>% pull(nm) %>% str_trim() %>% unique()
db_vocabs <- set_names(map(REPORT_DBS, vocab_of), REPORT_DBS)

audit <- imap_dfr(EUK_GROUPS, function(syn, g) {
  tibble(Group = g, Synonyms = length(syn),
         `In reference` = sum(syn %in% ref_vocab)) %>%
    bind_cols(as_tibble(map(db_vocabs, ~ sum(syn %in% .x))))
}) %>%
  mutate(Status = case_when(
    `In reference` == 0 ~ "absent from reference",
    if_any(all_of(REPORT_DBS), ~ .x == 0) ~ "MISSING IN >=1 DB - would read as empty",
    TRUE ~ "ok"
  ))

show_tbl(audit, "Mapping audit: how many of each group's declared synonyms exist in each vocabulary. A group matching nothing in a database silently reads as biologically absent - this is what happened to diatoms before `Diatomeae` was added. It also exposes an asymmetry that matters for the table below: the reference's vocabulary is NCBI/PR2-derived, so nr's names match it more readily than RepDB's UniEuk names do.")
```



Table: Mapping audit: how many of each group's declared synonyms exist in each vocabulary. A group matching nothing in a database silently reads as biologically absent - this is what happened to diatoms before `Diatomeae` was added. It also exposes an asymmetry that matters for the table below: the reference's vocabulary is NCBI/PR2-derived, so nr's names match it more readily than RepDB's UniEuk names do.

|Group                | Synonyms| In reference| clustnr| clustrepdb| nr| repdb|Status |
|:--------------------|--------:|------------:|-------:|----------:|--:|-----:|:------|
|Dinophyceae          |        3|            2|       1|          2|  1|     2|ok     |
|Ciliophora           |        1|            1|       1|          1|  1|     1|ok     |
|Apicomplexa+rel.     |        3|            3|       2|          2|  2|     2|ok     |
|Haptophyta           |        3|            3|       3|          2|  3|     2|ok     |
|Cryptophyta          |        3|            3|       2|          2|  2|     2|ok     |
|Bacillariophyta      |        5|            2|       4|          1|  3|     1|ok     |
|Pelagophyceae        |        1|            1|       1|          1|  1|     1|ok     |
|Dictyochophyceae     |        1|            1|       1|          1|  1|     1|ok     |
|other Stramenopiles  |        5|            2|       5|          4|  5|     4|ok     |
|Chlorophyta          |        6|            5|       4|          5|  4|     5|ok     |
|other Archaeplastida |        5|            4|       2|          3|  2|     3|ok     |
|Rhizaria             |        5|            4|       4|          4|  4|     4|ok     |
|Opisthokonta         |       11|           10|       9|          9|  9|     9|ok     |
|Amoebozoa            |        4|            4|       1|          2|  1|     2|ok     |
|Discoba/Excavata     |        5|            5|       4|          4|  4|     4|ok     |
|other Alveolata      |        1|            1|       1|          1|  1|     1|ok     |


``` r
euk_ref_grp <- ref %>%
  filter(str_starts(replace_na(ref_lineage, ""), "Eukaryota")) %>%
  transmute(seq, db = REF, grp = assign_euk_group(ref_lineage))

euk_db_grp <- lca %>%
  left_join(distinct(lin, db, taxid, lineage), by = c("db", "taxid")) %>%
  filter(str_detect(replace_na(lineage, ""), "Eukaryota")) %>%
  transmute(seq, db, grp = assign_euk_group(lineage))

euk_all <- bind_rows(euk_ref_grp, euk_db_grp)

# The share of eukaryotic calls the panel cannot place differs between sources
# (it depends on which intermediate nodes each taxonomy names), so the group
# percentages are renormalised over the assigned transcripts only - otherwise a
# source with more unplaced calls looks depleted in every group at once.
euk_unres <- euk_all %>%
  group_by(db) %>%
  summarise(`Eukaryotic calls` = n(),
            `Not placed by the panel` = sum(grp == UNRES),
            `%` = `Not placed by the panel` / `Eukaryotic calls` * 100,
            .groups = "drop")

euk_comp <- euk_all %>%
  filter(grp != UNRES) %>%
  count(db, grp) %>%
  group_by(db) %>% mutate(pct = n / sum(n) * 100) %>% ungroup()

N$n_euk_calls_repdb <- sum(euk_db_grp$db == "repdb")
N$n_euk_calls_ref <- nrow(euk_ref_grp)

show_tbl(euk_unres, "Eukaryotic calls the group panel cannot place. This is a property of the naming schemes, not of the biology, and is why the composition below is renormalised over placed transcripts.")
```



Table: Eukaryotic calls the group panel cannot place. This is a property of the naming schemes, not of the biology, and is why the composition below is renormalised over placed transcripts.

|db         | Eukaryotic calls| Not placed by the panel|    %|
|:----------|----------------:|-----------------------:|----:|
|Cohen2024  |           13,034|                   2,510| 19.3|
|clustnr    |            8,802|                   1,383| 15.7|
|clustrepdb |           11,170|                   2,694| 24.1|
|nr         |            8,886|                   1,424| 16.0|
|repdb      |           11,063|                   2,620| 23.7|

``` r
show_tbl(
  euk_comp %>% select(-n) %>%
    pivot_wider(names_from = db, values_from = pct) %>%
    arrange(desc(.data[[REF]])) %>%
    rename(Group = grp),
  sprintf("Eukaryotic community composition, as %% of the transcripts the panel places (%s eukaryotic transcripts in the reference, %s in RepDB).",
          fmt_n(N$n_euk_calls_ref), fmt_n(N$n_euk_calls_repdb))
)
```



Table: Eukaryotic community composition, as % of the transcripts the panel places (13,034 eukaryotic transcripts in the reference, 11,063 in RepDB).

|Group                | Cohen2024| clustnr| clustrepdb|   nr| repdb|
|:--------------------|---------:|-------:|----------:|----:|-----:|
|Dinophyceae          |      52.9|    49.5|       54.6| 48.9|  54.1|
|Haptophyta           |      13.8|    13.6|       13.0| 13.7|  13.0|
|Opisthokonta         |       6.1|     8.5|        3.5|  9.1|   3.4|
|Ciliophora           |       4.0|     2.5|        2.7|  2.4|   2.8|
|other Alveolata      |       3.5|     2.7|        4.5|  2.5|   5.3|
|other Stramenopiles  |       3.5|     5.2|        3.8|  5.3|   3.9|
|Chlorophyta          |       3.2|     4.0|        2.3|  4.0|   2.3|
|Pelagophyceae        |       2.5|     4.1|        2.7|  4.0|   2.7|
|Rhizaria             |       2.3|     1.7|        3.1|  1.7|   3.0|
|Bacillariophyta      |       2.2|     2.3|        1.6|  2.3|   1.5|
|Dictyochophyceae     |       1.5|     0.1|        1.1|  0.1|   1.0|
|Cryptophyta          |       1.4|     0.9|        1.0|  0.9|   1.0|
|other Archaeplastida |       1.3|     1.2|        0.4|  1.4|   0.3|
|Discoba/Excavata     |       0.8|     2.1|        5.0|  2.2|   5.0|
|Amoebozoa            |       0.6|     0.2|        0.6|  0.2|   0.6|
|Apicomplexa+rel.     |       0.4|     1.3|        0.2|  1.3|   0.3|


``` r
# how far each database's community profile sits from the reference's:
# total absolute deviation over the panel (0 = identical profile)
wide <- euk_comp %>% select(db, grp, pct) %>%
  pivot_wider(names_from = db, values_from = pct, values_fill = 0)
ref_pct <- wide[[REF]]

dev_tbl <- map_dfr(setdiff(names(wide), c("grp", REF)), function(d) {
  delta <- wide[[d]] - ref_pct
  tibble(db = d,
         `Total absolute deviation (pp)` = sum(abs(delta)),
         `Largest single deviation` = wide$grp[which.max(abs(delta))])
}) %>% arrange(`Total absolute deviation (pp)`)

N$dev_repdb <- dev_tbl$`Total absolute deviation (pp)`[dev_tbl$db == "repdb"]
N$dev_nr <- dev_tbl$`Total absolute deviation (pp)`[dev_tbl$db == "nr"]

show_tbl(dev_tbl, "Distance of each database's eukaryotic community profile from the reference's, summed over the panel (percentage points; lower is closer).")
```



Table: Distance of each database's eukaryotic community profile from the reference's, summed over the panel (percentage points; lower is closer).

|db         | Total absolute deviation (pp)|Largest single deviation |
|:----------|-----------------------------:|:------------------------|
|clustrepdb |                          16.3|Discoba/Excavata         |
|repdb      |                          16.7|Discoba/Excavata         |
|clustnr    |                          17.8|Dinophyceae              |
|nr         |                          19.2|Dinophyceae              |


``` r
ord <- euk_comp %>% filter(db == REF) %>% arrange(pct) %>% pull(grp)

plot_euk <- euk_comp %>%
  mutate(grp = factor(grp, levels = ord),
         db = factor(db, levels = DB_LEVELS)) %>%
  ggplot(aes(pct, grp, color = db)) +
  geom_line(aes(group = grp), color = "grey80", linewidth = .4) +
  geom_point(size = 2) +
  scale_db() +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  labs(x = "% of the transcripts the panel places", y = "") +
  thm

ggsave(fig("FigS3_euk_composition.pdf"), plot_euk, width = 7.2, height = 4.5)
plot_euk
```

<img src="/Users/gmutti/Desktop/projects/RepDBmaker/results/qc/comparison_files/figure-html/fig-euk-composition-1.png" alt="" width="691.2" />


``` r
# nested calls (RepDB says "other Alveolata", reference says "Dinophyceae") are
# coarser, not wrong; only genuinely different branches count as conflicts
nested <- function(a, b) {
  (!is.na(EUK_PARENT[a]) & EUK_PARENT[a] == b) | (!is.na(EUK_PARENT[b]) & EUK_PARENT[b] == a)
}

conc <- euk_db_grp %>%
  filter(db == "repdb") %>% select(seq, repdb = grp) %>%
  inner_join(select(euk_ref_grp, seq, ref = grp), by = "seq") %>%
  mutate(verdict = case_when(
    repdb == UNRES | ref == UNRES ~ "one side unresolved",
    repdb == ref ~ "same group",
    nested(repdb, ref) ~ "compatible (one coarser)",
    TRUE ~ "different group"
  ))

gp <- function(d, g) {
  v <- euk_comp$pct[euk_comp$db == d & euk_comp$grp == g]
  if (length(v) == 0) NA_real_ else v
}

conc_tbl <- conc %>% count(verdict) %>% mutate(`%` = n / sum(n) * 100)
scored <- conc %>% filter(verdict != "one side unresolved")
N$euk_concordance <- mean(scored$verdict != "different group") * 100

show_tbl(conc_tbl, "Per-transcript concordance on the curated panel, for transcripts both sources call eukaryotic.")
```



Table: Per-transcript concordance on the curated panel, for transcripts both sources call eukaryotic.

|verdict                  |     n|    %|
|:------------------------|-----:|----:|
|compatible (one coarser) |   477|  4.4|
|different group          |   725|  6.7|
|one side unresolved      | 3,568| 33.1|
|same group               | 6,019| 55.8|

``` r
show_tbl(
  scored %>% filter(verdict == "different group") %>%
    count(ref, repdb, sort = TRUE) %>% head(10),
  "Remaining genuine group-level disagreements (top 10)."
)
```



Table: Remaining genuine group-level disagreements (top 10).

|ref          |repdb               |  n|
|:------------|:-------------------|--:|
|Dinophyceae  |Discoba/Excavata    | 89|
|Opisthokonta |Discoba/Excavata    | 41|
|Dinophyceae  |other Stramenopiles | 38|
|Haptophyta   |Dinophyceae         | 33|
|Ciliophora   |Discoba/Excavata    | 27|
|Chlorophyta  |Dinophyceae         | 23|
|Ciliophora   |Dinophyceae         | 20|
|Haptophyta   |Discoba/Excavata    | 18|
|Ciliophora   |Rhizaria            | 17|
|Opisthokonta |Dinophyceae         | 16|

Across the 7,221 transcripts where both sources assign a group, they are consistent **90.0%** of the time. The composition profiles track each other closely: the reference's three dominant groups are also RepDB's, and no group is inverted.

One asymmetry to keep in mind: the panel is built from names, and the reference's vocabulary is NCBI/PR2-derived, so nr's lineages hit it more readily than RepDB's UniEuk ones - which is why a larger share of RepDB's eukaryotic calls goes unplaced. That is nomenclature, not biology, which is why the composition is renormalised over placed transcripts.

RepDB's community profile is also closer to the reference's than nr's overall (16.7 vs 19.2 percentage points of total deviation). Two divergences are worth stating plainly rather than averaging away: RepDB under-reports **Opisthokonta** (3.4% vs the reference's 6.1%), consistent with the downsampling cap measured in section 5, and over-reports **Discoba/Excavata** (5.0% vs 0.8%), which tracks the residual group-level conflicts above and looks like a UniEuk placement artifact rather than a coverage effect. nr, for its part, essentially misses **Dictyochophyceae** (0.1% vs 1.5%).

## 5. Is the gap concentrated in the clades RepDB caps?

This is the one place where the downsampling does show a cost. The caps in
`config/repdb.yaml` are Opisthokonta
(class, n=20), Ciliophora (order, n=20) and Embryophyta (family, n=20).


``` r
CAPPED <- c("Opisthokonta", "Ciliophora", "Embryophyta")

gap <- reason %>%
  filter(db == "repdb") %>%
  inner_join(ref, by = "seq") %>%
  filter(!is.na(ref_taxon), ref_taxon != "unclassified")

N$gap_total <- sum(gap$unresolved)
N$gap_pct <- mean(gap$unresolved) * 100

capped_tbl <- map_dfr(CAPPED, function(cl) {
  g <- gap %>% mutate(inside = str_detect(replace_na(ref_lineage, ""), cl))
  tibble(
    Clade = cl,
    `Queries in clade` = sum(g$inside),
    `% unresolved inside` = mean(g$unresolved[g$inside]) * 100,
    `% unresolved outside` = mean(g$unresolved[!g$inside]) * 100
  )
}) %>% mutate(Enrichment = `% unresolved inside` / `% unresolved outside`)

N$opis_in  <- capped_tbl$`% unresolved inside`[capped_tbl$Clade == "Opisthokonta"]
N$opis_out <- capped_tbl$`% unresolved outside`[capped_tbl$Clade == "Opisthokonta"]
N$capped_cost <- sum(capped_tbl$`Queries in clade` * capped_tbl$`% unresolved inside` / 100)

show_tbl(capped_tbl, "Are the transcripts RepDB fails to resolve concentrated in the downsampled clades?", digits = 2)
```



Table: Are the transcripts RepDB fails to resolve concentrated in the downsampled clades?

|Clade        | Queries in clade| % unresolved inside| % unresolved outside| Enrichment|
|:------------|----------------:|-------------------:|--------------------:|----------:|
|Opisthokonta |              642|               21.81|                13.52|       1.61|
|Ciliophora   |              420|               17.86|                13.70|       1.30|
|Embryophyta  |               11|               45.45|                13.77|       3.30|


``` r
show_tbl(
  gap %>% filter(unresolved) %>%
    count(`Reference rank` = ref_rank, sort = TRUE) %>%
    mutate(`%` = n / sum(n) * 100),
  "How deeply had the reference itself placed the transcripts RepDB leaves unresolved?"
)
```



Table: How deeply had the reference itself placed the transcripts RepDB leaves unresolved?

|Reference rank |   n|    %|
|:--------------|---:|----:|
|class          | 799| 29.7|
|order          | 605| 22.5|
|supergroup     | 446| 16.6|
|family         | 428| 15.9|
|division       | 189|  7.0|
|genus          | 169|  6.3|
|species        |  54|  2.0|

``` r
show_tbl(
  gap %>%
    mutate(ref_supergroup = replace_na(ref_supergroup, "(none)")) %>%
    count(ref_domain, ref_supergroup, unresolved) %>%
    pivot_wider(names_from = unresolved, values_from = n, values_fill = 0,
                names_prefix = "u") %>%
    transmute(Domain = ref_domain, Supergroup = ref_supergroup,
              Queries = uFALSE + uTRUE, Unresolved = uTRUE,
              `% unresolved` = uTRUE / (uFALSE + uTRUE) * 100) %>%
    arrange(desc(Unresolved)) %>% head(15),
  "Composition of the RepDB gap by reference supergroup (top 15)."
)
```



Table: Composition of the RepDB gap by reference supergroup (top 15).

|Domain    |Supergroup                   | Queries| Unresolved| % unresolved|
|:---------|:----------------------------|-------:|----------:|------------:|
|Eukaryota |Alveolata                    |   6,399|      1,023|         16.0|
|Eukaryota |(none)                       |   2,487|        374|         15.0|
|Bacteria  |Proteobacteria               |   3,552|        298|          8.4|
|Eukaryota |Hacrobia                     |   1,623|        278|         17.1|
|Eukaryota |Stramenopiles                |   1,020|        158|         15.5|
|Eukaryota |Opisthokonta                 |     642|        140|         21.8|
|Eukaryota |Archaeplastida               |     474|         73|         15.4|
|Bacteria  |(none)                       |     691|         62|          9.0|
|Bacteria  |Cyanobacteria                |     560|         42|          7.5|
|Bacteria  |Bacteroidetes/Chlorobi group |     399|         33|          8.3|
|Eukaryota |Rhizaria                     |     240|         33|         13.8|
|Bacteria  |Planctomycetes               |     171|         22|         12.9|
|Bacteria  |Actinobacteria               |     121|         19|         15.7|
|Bacteria  |Firmicutes                   |     124|         17|         13.7|
|Eukaryota |Excavata                     |      88|         16|         18.2|

Opisthokonta is measurably worse (21.8% unresolved inside versus 13.5% outside), and Ciliophora mildly so, which is consistent with the caps having a cost. But the affected clades are a small part of a marine metatranscriptome: the three capped clades together account for about 220 unresolved transcripts, ~1.1% of the query set. And the transcripts RepDB fails to place are overwhelmingly ones the reference itself placed only at class or order level - not well-characterised sequences being lost.

# Session info


``` r
sessionInfo()
```

```
## R version 4.5.2 (2025-10-31)
## Platform: aarch64-apple-darwin20
## Running under: macOS Sequoia 15.6
## 
## Matrix products: default
## BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
## 
## locale:
## [1] C/UTF-8/C/C/C/C
## 
## time zone: Europe/Madrid
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] patchwork_1.3.2 lubridate_1.9.4 forcats_1.0.1   stringr_1.6.0  
##  [5] dplyr_1.1.4     purrr_1.2.2     readr_2.1.6     tidyr_1.3.2    
##  [9] tibble_3.3.0    ggplot2_4.0.1   tidyverse_2.0.0
## 
## loaded via a namespace (and not attached):
##  [1] sass_0.4.10        generics_0.1.4     stringi_1.8.7      hms_1.1.4         
##  [5] digest_0.6.39      magrittr_2.0.4     evaluate_1.0.5     grid_4.5.2        
##  [9] timechange_0.3.0   RColorBrewer_1.1-3 bookdown_0.46      fastmap_1.2.0     
## [13] jsonlite_2.0.0     scales_1.4.0       textshaping_1.0.4  jquerylib_0.1.4   
## [17] cli_3.6.5          rlang_1.3.0        crayon_1.5.3       bit64_4.6.0-1     
## [21] withr_3.0.2        cachem_1.1.0       yaml_2.3.12        otel_0.2.0        
## [25] tools_4.5.2        parallel_4.5.2     tzdb_0.5.0         vctrs_0.7.3       
## [29] R6_2.6.1           lifecycle_1.0.4    bit_4.6.0          vroom_1.6.7       
## [33] ragg_1.5.0         pkgconfig_2.0.3    pillar_1.11.1      bslib_0.9.0       
## [37] gtable_0.3.6       glue_1.8.0         rmdformats_1.0.4   systemfonts_1.3.1 
## [41] xfun_0.55          tidyselect_1.2.1   knitr_1.51         farver_2.1.2      
## [45] htmltools_0.5.9    labeling_0.4.3     rmarkdown_2.30     compiler_4.5.2    
## [49] S7_0.2.1
```
