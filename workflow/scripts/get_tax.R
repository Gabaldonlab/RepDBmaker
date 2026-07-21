# Harmonize a source taxonomy (UniProt NCBI lineage or P10K native lineage) into
# the UniEuk / EukProt framework.
#
# Idea (see test/tax_harmony/test.R): EukProt defines, for its genera, which
# UniEuk clade names sit at each output rank:
#     sk (Taxonomy_UniEuk[2], e.g. Amorphea) -> p     Supergroup_UniEuk -> c
#     Taxogroup1_UniEuk                       -> o     Taxogroup2_UniEuk -> f
# Using those vocabularies we classify every UniEuk lineage into p/c/o/f, keyed
# by the clade it ends at. Then, for each proteome, we look up its genus in that
# table; if the genus is not in UniEuk we walk one rank up (family, order, ...)
# and look up that clade instead. The genus/species themselves come from the
# proteome. This is a vectorized lookup + join (no per-proteome regex scan).
suppressMessages(library(tidyverse))

# ---- inputs ------------------------------------------------------------------
lineage <- read_delim(snakemake@input[["lineage"]],
  col_names = FALSE, show_col_types = FALSE, delim = "\t",
  col_types = cols(.default = "c")  # taxid stays character (UniProt) for if_else
)
unieuk <- read_delim(snakemake@input[["unieuk"]], show_col_types = FALSE)
ep <- read_delim(snakemake@input[["ep"]], show_col_types = FALSE) %>%
  select(-Replaces_EukProt_ID)
notep <- read_delim(snakemake@input[["notep"]], show_col_types = FALSE) %>%
  select(-Replaced_by_EukProt_ID)

# UniEuk wraps many clade names in single quotes (e.g. 'coccidiomorphea', 'core
# centramoebids') while EukProt's vocabularies do not, so names must be
# normalized (unquoted, trimmed) before any comparison, or valid families are
# missed (coccidiomorphea).
norm <- function(x) trimws(gsub("'", "", x))

# ---- EukProt vocabularies: UniEuk clade name -> output rank -------------------
ep_all <- bind_rows(ep, notep) %>%
  mutate(sk = map_chr(str_split(Taxonomy_UniEuk, ";"), 2))
sk_vocab  <- na.omit(unique(norm(ep_all$sk)))                 # -> p (supergroup, e.g. Amorphea)
sup_vocab <- na.omit(unique(norm(ep_all$Supergroup_UniEuk)))  # -> c
tg1_vocab <- na.omit(unique(norm(ep_all$Taxogroup1_UniEuk)))  # -> o
tg2_vocab <- na.omit(unique(norm(ep_all$Taxogroup2_UniEuk)))  # -> f

# a clade name -> its output rank slot (first match wins: p > c > o > f)
clade_rank <- c(
  setNames(rep("p", length(sk_vocab)),  sk_vocab),
  setNames(rep("c", length(sup_vocab)), sup_vocab),
  setNames(rep("o", length(tg1_vocab)), tg1_vocab),
  setNames(rep("f", length(tg2_vocab)), tg2_vocab)
)
clade_rank <- clade_rank[!duplicated(names(clade_rank))]

# ---- classify every UniEuk lineage into p/c/o/f, keyed by its last clade ------
clades_raw <- lapply(strsplit(unieuk$taxon, ";"), function(x) {
  x <- norm(x)
  x[!is.na(x) & x != ""]
})

# the clade each UniEuk entry ends at (its key)
last_vec <- vapply(clades_raw, function(x) if (length(x)) x[length(x)] else NA_character_, character(1))
last_clade <- tibble(index = seq_along(last_vec), clade = last_vec) %>% filter(!is.na(clade))

# long form: one row per (entry, level); the deepest level per rank slot wins
uni_long <- tibble(
  index = rep.int(seq_along(clades_raw), lengths(clades_raw)),
  pos   = sequence(lengths(clades_raw)),
  value = unlist(clades_raw, use.names = FALSE),
  slot  = unname(clade_rank[unlist(clades_raw, use.names = FALSE)])
) %>%
  filter(!is.na(slot))

slots <- uni_long %>%
  arrange(index, desc(pos)) %>%
  distinct(index, slot, .keep_all = TRUE) %>%   # first = deepest level for that slot
  pivot_wider(id_cols = index, names_from = slot, values_from = value)

unieuk_clade_lineage <- last_clade %>%
  left_join(slots, by = "index") %>%
  select(clade, any_of(c("p", "c", "o", "f")))
for (col in c("p", "c", "o", "f")) {
  if (!col %in% names(unieuk_clade_lineage)) unieuk_clade_lineage[[col]] <- NA_character_
}
# drop clades whose lineage is ambiguous (same name under conflicting
# supergroups, e.g. the animal "Vertebrata" vs the red-alga genus "Vertebrata")
unieuk_clade_lineage <- unieuk_clade_lineage %>%
  distinct() %>%
  group_by(clade) %>% filter(n() == 1) %>% ungroup()

# EukProt's explicit genus records: richer than the classified UniEuk path
# because EukProt assigns Taxogroup1/Taxogroup2 (order/family) that the UniEuk
# lineage does not always spell out (e.g. Plasmodium -> coccidiomorphea).
eukprot_genus_lineage <- ep_all %>%
  transmute(
    clade = norm(gsub("_", "", Genus_UniEuk)),
    p = norm(sk), c = norm(Supergroup_UniEuk),
    o = norm(Taxogroup1_UniEuk), f = norm(Taxogroup2_UniEuk)
  ) %>%
  filter(!is.na(clade), clade != "") %>%
  distinct(clade, .keep_all = TRUE)

manual <- tribble(
  ~clade,          ~p,               ~c,             ~o,             ~f,
  "Oomycota",      "Diaphoretickes", "Stramenopiles", "other_Gyrista", "Peronosporomycetes",
  "Heterolobosea", "Discoba",        "Heterolobosea", "Heterolobosea", "Heterolobosea",
  # Chordata in the (NCBI) lineage is the unambiguous animal signal, so resolve
  # the "Vertebrata" homonym (animal subphylum vs red-alga genus) in its favour.
  # The non-vertebrate chordate subphyla are more specific candidates (they sit
  # below Chordata in the lineage) so they win the family over "Vertebrata".
  "Chordata",         "Amorphea", "Opisthokonta", "Metazoa", "Vertebrata",
  "Cephalochordata",  "Amorphea", "Opisthokonta", "Metazoa", "Cephalochordata",
  "Tunicata",         "Amorphea", "Opisthokonta", "Metazoa", "Urochordata",
  "Urochordata",      "Amorphea", "Opisthokonta", "Metazoa", "Urochordata"
)

# priority when a clade name appears in several sources: manual > EukProt genus
# record > UniEuk-classified clade
clade_lineage <- bind_rows(manual, eukprot_genus_lineage, unieuk_clade_lineage) %>%
  distinct(clade, .keep_all = TRUE)

# ---- per-proteome candidate taxa (most specific -> most general) -------------
# clean each token: strip the P10K rank prefix ("P_", ...) and any "_X" incertae
# sedis suffix PER TOKEN (cleaning the whole string truncated lineages at an
# internal underscore, dropping the genus). Apply the small manual renames.
clean_tokens <- function(s) {
  t <- strsplit(s, ";")[[1]]
  t <- norm(gsub("_.*", "", gsub("^[A-Z]_", "", t)))
  t <- recode(t,
    Fragilariophyceae = "Diatomeae", Thalassiosiraceae = "Diatomeae",
    Echinamoebida = "Echinamoebidia", Haemosporida = "Haemospororida"
  )
  t[!is.na(t) & !t %in% c("cellular organisms", "Eukaryota", "root", "")]
}
# the species name (last lineage element); keep it as-is apart from a stray rank
# prefix - do NOT strip "_..." here or strain IDs get mangled (ST99CH_3D7).
last_token <- function(s) {
  t <- strsplit(s, ";")[[1]]
  gsub("^[A-Z]_", "", t[length(t)])
}

prot <- lineage %>%
  transmute(mnemo = X1, taxid = X2, raw = X3) %>%
  mutate(
    # P10K carries the species in column 2 (non-numeric); UniProt a numeric taxid
    species = if_else(grepl("\\D", taxid), taxid, map_chr(raw, last_token)),
    genus   = word(species, 1),
    cand    = map(raw, ~ rev(clean_tokens(.x)))  # specific -> general
  )

# For each proteome, join every candidate to the lookup and fill each rank from
# the most specific candidate that actually knows it. (A single most-specific
# match is not enough: a species-level UniEuk entry may match but only carry the
# supergroup, while the genus entry carries the order/family.)
matched <- prot %>%
  select(mnemo, cand) %>%
  unnest_longer(cand, indices_to = "priority") %>%
  inner_join(clade_lineage, by = c("cand" = "clade")) %>%
  mutate(across(c(p, c, o, f), ~ na_if(., ""))) %>%
  arrange(mnemo, priority) %>%   # most specific candidate first
  group_by(mnemo) %>%
  summarise(across(c(p, c, o, f), ~ first(na.omit(.))), .groups = "drop")

result <- prot %>%
  select(mnemo, genus, species) %>%
  inner_join(matched, by = "mnemo")

# fill down: an empty rank inherits the nearest known ancestor, marked
# "unassigned_" (so every lineage is a complete 7 ranks and the taxdump gets
# meaningful nodes instead of blanks). p/c/o/f only; genus/species always set.
rk <- c("p", "c", "o", "f")
m <- as.matrix(result[rk])
m[m == ""] <- NA
last_real <- rep("Eukaryota", nrow(m)) # base ancestor = the domain
for (j in seq_along(rk)) {
  empty <- is.na(m[, j])
  m[empty, j] <- paste0("unassigned_", last_real[empty])
  last_real[!empty] <- m[!empty, j]
}
for (j in seq_along(rk)) result[[rk[j]]] <- m[, j]

unmatched <- setdiff(prot$mnemo, matched$mnemo)
if (length(unmatched)) {
  message(length(unmatched), " proteomes had no UniEuk clade in their lineage.")
}

# ---- write d__..;s__.. lineages ----------------------------------------------
blank <- function(x) coalesce(x, "")
new_taxonomy <- result %>%
  transmute(
    sp = mnemo,
    d = "d__Eukaryota",
    p = paste0("p__", blank(p)), c = paste0("c__", blank(c)),
    o = paste0("o__", blank(o)), f = paste0("f__", blank(f)),
    g = paste0("g__", blank(genus)), s = paste0("s__", blank(species))
  ) %>%
  arrange(p, c, f, o, g) %>%
  unite(taxon, d, p, c, o, f, g, s, sep = ";")

write.table(new_taxonomy,
  quote = FALSE, sep = "\t", row.names = FALSE,
  col.names = FALSE, snakemake@output[[1]]
)
