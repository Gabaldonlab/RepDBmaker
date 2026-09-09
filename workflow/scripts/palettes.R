# Every palette below that reaches a paper figure is verified by
# workflow/scripts/check_palettes.R: OKLCH lightness band, chroma floor, WCAG
# contrast on white, and OKLab Delta E between all pairs under normal,
# protanope and deuteranope vision (Machado-Oliveira-Fernandes 2009, severity
# 1.0). Run that script after touching anything here - do not eyeball it.

# Kingdoms. Purple sits at L=0.44 rather than mid-lightness: at L=0.55 it
# collapsed against Archaea blue (dE 7.0 deuteranope, below the 8 floor).
color_kingdoms <- c("#D9043D", "#56B4E9", "#009E73", "#CC79A7")
names(color_kingdoms) <- c("Eukaryota", "Archaea", "Bacteria", "Viruses")

# "Ambig" is an absence category, so it is deliberately neutral and is the one
# slot exempt from the chroma floor; it still clears dE 15 against the rest.
color_contaminants <- c("#B4B4B4", "#418AD1", "#50A064", "#70279A")
names(color_contaminants) <- c("Ambig", "Archaea", "Bacteria", "Viruses")

# Benchmarked databases: hue carries the family (nr blue, RepDB green) and
# lightness carries clustered-or-not, so the pairing is readable without a
# fifth and sixth hue. Both members of each pair sit inside the L 0.43-0.77
# band - the previous "clustrepdb" (#62E4AB, L=0.83) sat above it and rendered
# at 1.55:1 on white, i.e. invisible when printed.
color_benchdbs <- c("nr"="#285994",
                    "clustnr"="#689BDB",
                    "repdb"="#267B4C",
                    "clustrepdb"="#6EBF8C")

# Panel C of Fig. 5: agreement is green, disagreement orange, and the two
# "no answer" outcomes are neutral - light grey where a homolog exists but the
# LCA would not commit, dark grey where there is no homolog at all. Both greys
# are absence categories and so waive the chroma floor; they still separate
# from each other by dE 31 (they differ by 0.31 in OKLCH L).
#
# The green is #349D62 rather than a deeper forest green because it also has to
# stay clear of `repdb`'s green in panels A/B/D of the same figure: at
# L=0.55/H=150 the two were dE 3.6 apart, which is indistinguishable. Pushing it
# lighter opens that to dE 10.2 while keeping green-vs-orange above the dE 8 CVD
# floor (8.7) - a yellower green scores better against repdb but collapses
# against the orange under protanopia (dE 0.4).
outcome_colors_ref <- c("#349D62", "#CF6F19", "#B4B4B4", "#585858")

color_dbs <- c("p10k"="#449E77",
               "eukprot"="#EBCD62",
               "uniprot"="#CA3D54",
               "custom"="#933E9A")

# full source palette: the four eukaryotic dbs plus prokaryotes/viruses/other
color_sources <- c(color_dbs, "gtdb"="#3d6dcc", "virus"="#7a7a7a", "other"="#cccccc")

type_palette <- c("Genome"="#762A83",
                  "S.C. genome"="#C2A5CF",
                  "Transcriptome"="#1B7837",
                  "S.C. transcriptome"="#D9F0D3",
                  # "ranscriptome,EST"="#5AAE61",
                  "EST"="#5AAE61")


ass_palette <- c(
  "Complete genome" = "#08306B", # Darkest Navy
  "Chromosome"      = "#2171B5", # Strong Royal Blue
  "Scaffold"        = "#6BAED6", # Mid Sky Blue
  "Contig"          = "#C6DBEF" # Pale Blue
  # "NA"              = "grey80"
)


busco_colors <- c("#2CBBEF","#0099CF","#F3E600","#FF343E")
names(busco_colors) <- c("Single", "Duplicated", "Fragmented", "Missing")

busco_colors2 <- c("#2CBBEF","#0099CF","#F3E600","#FF343E")
names(busco_colors2) <- c("Complete", "Duplicated", "Fragmented", "Missing")


phylum_col <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999")
names(phylum_col) <- c("Amorphea", "Ancyromonadida","CRuMs","Diaphoretickes","Discoba","Malawimonadida","Metamonada")

omark_complete_cols <- c("#46bea1ff","#cfee8eff","#ed1c5aff")
names(omark_complete_cols) <- c("Single", "Duplicated","Missing")

omark_consistent_cols <- c("#000000ff", "#8a5df1ff", "#e69f00ff", "#3db7e9ff")
names(omark_consistent_cols) <- c("Unknown","Inconsistent","Contaminant","Consistent")

tiara_palette <- c(
  "grey", "#C3461F", "#437848",
  "#2DA0B0", "#92CAC1", "#114E76",
  "grey40", "black"
)
names(tiara_palette) <- c(
  "eukarya", "mitochondrion", "plastid",
  "archaea", "bacteria", "prokarya",
  "short", "unknown"
)

omark_kingdom_palette <- c("Eukaryota"="grey",
                           "Bacteria"="#92CAC1", 
                           "Archaea"="#114E76",
                           "Ambiguous contaminant"="black")


euka_groups <- c("#F2BF40", "#8AC4DB", "#F68512", "#BB3F3B",
                        "#6FBE47", "#B38481", "#908DCD",
                        "#8F97A0", "grey80", "grey80", "grey80",
                        "#A8995C", "grey80", "#865F84")
names(euka_groups) <- c("Stramenopiles", "Rhizaria", "Alveolata", "Telonemia",
                 "Archaeplastida", "Haptista", "Cryptista",
                 "Amorphea", "Discoba", "Malawimonadida", "Metamonada",
                 "CRuMs", "Ancyromonadida", "Hemimastigophora")

colors_euka_phylum <- c("Diaphoretickes" = "#F2BF40",
                        "Amorphea" = "#8AC4DB", 
                        "Metamonada" = "#F68512", 
                        "Discoba" = "#BB3F3B",
                        "CRuMs" = "#6FBE47", 
                        "Ancyromonadida" = "#B38481", 
                        "Malawimonadida" = "#908DCD")

virus_colors <- c(
  # --- DNA Groups (Class I & II) ---
  "DNA"          = "#50728B", # Muted blue (Generic)
  "dsDNA"        = "#2A5676", # Dark navy (Class I)
  "ssDNA"        = "#709AE1", # Bright blue (Class II)
  "ssDNA(+)"     = "#99C1F1", # Light blue
  "ssDNA(-)"     = "#4A90E2", # Mid blue
  "ssDNA(+/-)"   = "#5C7EB5", # Steel blue
  "dsDNA; ssDNA" = "#63B3ED", # Cyan-blue mix
  
  # --- RNA Groups (Class III, IV & V) ---
  "RNA"          = "#B56576", # Dusty rose (Generic)
  "dsRNA"        = "#6D213C", # Deep burgundy (Class III)
  "ssRNA"        = "#F28123", # Orange (Generic ssRNA)
  "ssRNA(+)"     = "#F29F05", # Golden orange (Class IV)
  "ssRNA(-)"     = "#D9525E", # Coral red (Class V)
  "ssRNA(+/-)"   = "#E57C51", # Burnt orange
  
  # --- Reverse Transcribing (Class VI & VII) ---
  "ssRNA-RT"     = "#4C956C", # Leaf green (Class VI)
  "dsDNA-RT"     = "#2C6E49", # Dark forest green (Class VII)
  
  # --- Unknowns ---
  "unknown"      = "#A0A0A0"
  # "NA"           = "#D3D3D3"
)

gff_columns <- c("seqid", "source", "type", "start", "end", 
                 "score", "strand", "phase", "attributes")


theme_legend <- theme(axis.text.x = element_text(angle=45, size=7, hjust=1),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.subtitle = element_text(vjust=-.9, size=9))
