# Load required libraries
suppressMessages(library(tidyverse))  # Load the tidyverse package for data manipulation
# library(optparse)  # Load the optparse package for command-line argument parsing

# Read the p10k taxonomy data (see README)
p10k <- read_delim(snakemake@input[["p10k"]], col_names = T, 
                   show_col_types = FALSE, delim = "\t") %>% 
                   mutate(lineage=gsub("_X*;", ";", lineage))
# Read the UniEuk data from a fixed file
unieuk <- read_delim(snakemake@input[["unieuk"]], show_col_types = FALSE)

# Split the 'taxon' column of UniEuk data into separate taxonomic levels
clades_raw <- strsplit(unieuk$taxon, ";")

# Create a tibble (data frame) from the taxonomic levels, padding shorter sequences with NA values
taxonomy <- as_tibble(do.call(rbind, lapply(clades_raw, `length<-`, max(lengths(clades_raw)))))
# Extract the last taxonomic clade for each species
last_clade <- apply(taxonomy, 1, function(x) tail(na.omit(x), 1))

# Read EukProt data from multiple files and select relevant columns
eukprot_groups_included <- read_delim(snakemake@input[["ep"]], show_col_types = FALSE) %>%
    select(Supergroup_UniEuk, Taxogroup1_UniEuk, Taxogroup2_UniEuk) %>%
    distinct()

eukprot_groups_excluded <- read_delim(snakemake@input[["notep"]], show_col_types = FALSE) %>%
    select(Supergroup_UniEuk, Taxogroup1_UniEuk, Taxogroup2_UniEuk) %>%
    distinct()

eukprot_groups <- rbind(eukprot_groups_included, eukprot_groups_excluded)

# Initialize a new taxonomy table with specific columns
new_taxonomy <- tibble(sp = as.character(), id = as.numeric(), d = "Eukaryota", p = as.character(),
                       c = as.character(), o = as.character(), f = as.character(), 
                       g = as.character(), s = as.character())

# Loop through each species in the p10k taxonomy data
for (species in 1:nrow(p10k)) {
    mnemo <- pull(p10k[species, "p10k_id"])  # Get the mnemonic identifier
    # taxid <- pull(p10k[species, "X2"])  # Get the taxonomic identifier

    els <- str_split(gsub("[A-Z]_", "", p10k[species, ]$lineage), pattern = ";", simplify = T)
    last_tax <- els[length(els)]  # Get the last taxonomic level
    els <- els[2:length(els)]  # Exclude the first two levels
    # Split the taxonomic path

    found <- FALSE

    # Initialize variables for new taxonomic levels
    new_class <- ""
    new_order <- ""
    new_family <- ""
    new_sk <- ""
    # Iterate through taxonomic levels to find matches in UniEuk data
    for (taxon in rev(els)) {
        # Check if the taxon matches a clade in the UniEuk data
        sub_df <- taxonomy[which(grepl(paste0("^", taxon, "$"), last_clade)),]
        if (nrow(sub_df) > 0) {
            found <- TRUE
            # Extract taxonomic information from UniEuk data
            new_class <- sub_df[which(sub_df %in% eukprot_groups$Supergroup_UniEuk)]
            new_class <- ifelse(ncol(new_class) > 0, pull(new_class), "")
            new_order <- sub_df[which(sub_df %in% eukprot_groups$Taxogroup1_UniEuk)]
            new_order <- ifelse(ncol(new_order) > 0, pull(new_order), "")
            new_family <- sub_df[which(sub_df %in% eukprot_groups$Taxogroup2_UniEuk)]
            new_family <- ifelse(ncol(new_family) > 0, pull(new_family), "")

            new_sk <- pull(sub_df[, 'V2'])

            # Add a new row to the taxonomy table with taxonomic information
            new_taxonomy <- add_row(new_taxonomy, sp = mnemo, 
            # id = taxid,
                                    d = paste0("d__Eukaryota"), p = paste0("p__",new_sk), 
                                    c = paste0("c__",new_class), o = paste0("o__", new_order),
                                    f = paste0("f__",new_family), g = paste0("g__",str_split(last_tax, " ")[[1]][1]), 
                                    s = paste0("s__", p10k[species,]$species))
            break
        }
    }
}

# Arrange the new taxonomy table by taxonomic levels and write it to the output file
new_taxonomy <- arrange(new_taxonomy, p, c, f, o, g) %>% 
    unite(taxon,d,p,c,o,f,g,s, sep = ";")  %>% 
    select(-id)

# write_delim(new_taxonomy, stdout(), delim = "\t", col_names = FALSE)
write.table(new_taxonomy, quote=FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, snakemake@output[[1]])