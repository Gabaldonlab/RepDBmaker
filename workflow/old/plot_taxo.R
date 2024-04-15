library(tidyverse)
library(data.tree)
library(networkD3)
library(ggraph)
library(igraph)

theme_set(theme_bw())

euka_groups <- c("Amorphea","Ancyromonadida","CRuMs","Diaphoretickes",
                 "Discoba","Hemimastigophora","Malawimonadida","Metamonada","Provora")

color_df <- tibble(group = euka_groups,
                   color = c("#F2BF40", "#8AC4DB", "#F68512", "#BB3F3B",
                             "#6FBE47", "#B38481", "#908DCD","#A8995C", "#865F84"))

taxo <- read_delim("data/taxonomies/uniprot_eukprot_taxonomy.tsv", 
                   col_names = c("mnemo","taxon"), delim = "\t") %>% 
  tidyr::separate(taxon, sep = ";", into = c("k", "p", "c", "o", "f", "g", "s")) %>% 
  mutate_all(~gsub("*.__", "", .x))

taxo_red <- taxo %>% 
  select(-c(mnemo,g,s)) %>% 
  distinct()

taxo_red$pathString <- paste(taxo_red$k, taxo_red$p, taxo_red$c, taxo_red$o, taxo_red$f, sep = "/")

taxo_network <- ToDataFrameNetwork(as.Node(taxo_red), "name") %>% 
  left_join(color_df, by=c("name"="group")) %>% 
  mutate(from = sub(".*\\/", "", from),
         to = sub(".*\\/", "", to),
         color = ifelse(is.na(color), "grey30", color)) %>% 
  distinct() %>% 
  filter(from != to)

# ggraph(taxo_network, 'partition', circular = TRUE) + 
#   geom_node_arc_bar(aes(fill = depth), size = 0.25) + 
#   coord_fixed() +
#   theme(legend.position = "none")

terminal <- names(which(degree(graph_from_data_frame(taxo_network), mode="out")<=0))

v_meta <- distinct(taxo_network[, c("name", "color")]) %>% 
  add_row(name="Eukaryota", color="red") %>% 
  rowwise() %>% 
  mutate(size = sum(apply(taxo, 1, function(r) any(r %in% name))),
         new_lab = ifelse(name %in% color_df$group, paste(name, "-",size), NA),
         lab_terminal = ifelse(name %in% terminal & is.na(new_lab), name, NA))

g <- igraph::graph_from_data_frame(taxo_network, vertices = v_meta)

p1 <- ggraph(g, 'tree') + 
  geom_edge_diagonal() +
  geom_node_point(aes(fill=color, size=size), pch=21) +
  scale_y_reverse() +
  coord_flip(ylim = c(4,-.7)) +
  geom_node_label(aes(label = new_lab, fill=color), label.size = NA, color="white", family="Helvetica", 
                  repel = T, direction="both", force=3, fontface="bold") +
  geom_node_text(aes(label = lab_terminal), family="Helvetica", 
                 hjust=0, nudge_y=.05, fontface="bold", color="grey50") +
  scale_color_identity() +
  scale_fill_identity() +
  theme(legend.position = "none")


ggsave("plots/broaddb_family_euka.png", p1, dpi = 300, height = 12, width = 10)

# TODO
# maybe plot also complete unieuk taxonomy to see what changes!!!!