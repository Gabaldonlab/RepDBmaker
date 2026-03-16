color_kingdoms <- c("#D9043D", "#56B4E9", "#009E73", "#CC79A7")
names(color_kingdoms) <- c("Eukaryota", "Archaea", "Bacteria", "Viruses")

color_contaminants <- c("grey50", "#56B4E9", "#009E73", "#CC79A7")
names(color_contaminants) <- c("Ambig", "Archaea", "Bacteria", "Viruses")

color_benchdbs <- c("nr"="#20558A",
                    "clustnr"="#307FCF",
                    "repdb"="#449E77",
                    "clustrepdb"="#62E4AB")

color_dbs <- c("p10k"="#449E77", 
               "eukprot"="#EBCD62", 
               "uniprot"="#CA3D54",
               "custom"="#933E9A")

type_palette <- c("genome"="#762A83",
                  "single-cell genome"="#C2A5CF",
                  "transcriptome"="#1B7837",
                  "single-cell transcriptome"="#D9F0D3",
                  "transcriptome,EST"="#5AAE61",
                  "EST"="#5AAE61")

ass_palette <- c(
  "Complete Genome" = "#08306B", # Darkest Navy
  "Chromosome"      = "#2171B5", # Strong Royal Blue
  "Scaffold"        = "#6BAED6", # Mid Sky Blue
  "Contig"          = "#C6DBEF", # Pale Blue
  "NA"              = "grey80"
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
