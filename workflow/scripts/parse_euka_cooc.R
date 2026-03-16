tax <- read_delim("results/taxonomies/repdb.tsv",
                  col_names = c("mnemo", "k", "p", "c", "o", "f", "g", "s")
)

pair_cooc <- read_delim("results/contamination/pair_counts.tsv",
                        col_names = c("n", "tmpq", "tmpt")
) %>%
  mutate(
    q = pmin(tmpq, tmpt),
    t = pmax(tmpq, tmpt)
  ) %>%
  select(-contains("tmp")) %>%
  group_by(q, t) %>%
  summarise(n = sum(n))

euka_mat <- pair_cooc %>%
  # filter(n>1000, q!=t) %>%
  left_join(tax, by = c("q" = "mnemo")) %>%
  left_join(tax, by = c("t" = "mnemo")) %>%
  filter(k.x == "Eukaryota", k.y == "Eukaryota", p.x != p.y)
