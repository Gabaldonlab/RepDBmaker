repdb_counts <- repdb_df %>% 
  group_by(k, label) %>% 
  summarise(n = sum(num_seqs),
            # mn=mean(num_seqs), var=var(num_seqs), 
            .groups = "drop")

tm <- voronoiTreemap(
  data = filter(repdb_counts, n>100000),
  levels = c("k", "label"),
  cell_size = "n",
  shape = "rectangle",
  seed = 2
)

pastel_teal_vec <- c("#4dbf9f", "#80cdc1", "#b2dfdb", "#a1d99b")
pdf(NULL)
drawTreemap(tm, 
            label_level = NULL,
            color_palette = pastel_teal_vec,
            # color_palette =c("grey99", "grey70", "grey85", "grey90"),
            # color_palette = color_kingdoms[unique(test$k)],
            border_level = c(1,2),
            border_color = c("white", "white"),
            border_size = c(3, 1),
            width = .9, height = .9) 
p_tm_logo <- grid::grid.grab()
dev.off()

ggsave("test/logo/repdb_logo_exp.pdf", p_tm_logo, width = 6, height = 6*9/16)
ggsave("test/logo/repdb_logo_exp.pdf", p_tm_logo, width = 4, height = 4)
ggsave("test/logo/repdb_logo_exp.png", p_tm_logo, dpi = 600, width = 4, height = 4)
