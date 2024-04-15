# Load required libraries
suppressMessages(library(tidyverse))  # Load the tidyverse package for data manipulation

# Read the Original taxonomy data (see README)
ncbi <- read_delim(snakemake@input[["lineage"]], col_names = F, 
                   show_col_types = FALSE, delim = "\t")

# Read the UniEuk taxonomy
unieuk <- read_delim(snakemake@input[["unieuk"]], show_col_types = FALSE)

# Split the 'taxon' column of UniEuk data into separate taxonomic levels
clades_raw <- strsplit(unieuk$taxon, ";")

# Create a data frame from the taxonomic levels, padding shorter sequences with NA values
taxonomy <- as_tibble(do.call(rbind, lapply(clades_raw, `length<-`, max(lengths(clades_raw)))))
# Extract the last taxonomic clade for each species
last_clade <- apply(taxonomy, 1, function(x) tail(na.omit(x), 1))

# Read EukProt
ep <- read_delim(snakemake@input[["ep"]], show_col_types = FALSE) %>% 
    select(-Replaces_EukProt_ID)
notep <- read_delim(snakemake@input[["notep"]], show_col_types = FALSE) %>% 
    select(-Replaced_by_EukProt_ID)
ep <- rbind(ep, notep) %>% 
    rowwise() %>% 
    mutate(sk=str_split(Taxonomy_UniEuk, ";")[[1]][2])

eukprot_groups <- ep %>%
    select(sk, Supergroup_UniEuk, Taxogroup1_UniEuk, Taxogroup2_UniEuk, Genus_UniEuk) %>%
    distinct()

# Initialize a new taxonomy table with specific columns
new_taxonomy <- tibble(sp = as.character(), d = "Eukaryota", p = as.character(),
                       c = as.character(), o = as.character(), f = as.character(), 
                       g = as.character(), s = as.character())

# Loop through each species in the NCBI taxonomy data
for (species in 1:nrow(ncbi)) {
    mnemo <- pull(ncbi[species, "X1"])  # Get the mnemonic identifier
    taxid <- pull(ncbi[species, "X2"])  # Get the taxonomic identifier
    els <- str_split(gsub("_.*", "", gsub("[A-Z]_", "", pull(ncbi[species, 3]))), 
                     pattern = ";", simplify = T) # Split the taxonomic path
    last_tax <- els[length(els)]  # Get the last taxonomic level
    els <- els[3:length(els)]  # Exclude the first two levels
    
    if (grepl("\\D", taxid)) {
        species_name <- taxid # this means that it comes from p10k as second column is species
    } else {
        species_name <- last_tax
    }
    genus <- str_split(species_name, " ")[[1]][1]
    # If the genus is already in Eukprot then it's easy
    if (genus %in% gsub("_", "", ep$Genus_UniEuk)){
        groups <- as.character(eukprot_groups[eukprot_groups$Genus_UniEuk==genus,])
        new_taxonomy <- add_row(new_taxonomy, sp = mnemo, 
                                d = paste0("d__Eukaryota"), p = paste0("p__", groups[1]), 
                                c = paste0("c__", groups[2]), o = paste0("o__", groups[3]),
                                f = paste0("f__", groups[4]), g = paste0("g__", genus), 
                                s = paste0("s__", species_name))
    } else {
        # Fragilariophyceae and Thalassiosiraceae to Diatomeae
        els <- gsub("Fragilariophyceae", "Diatomeae", els)
        els <- gsub("Thalassiosiraceae", "Diatomeae", els)
        # And some other small tweaks
        els <- gsub("Echinamoebida", "Echinamoebidia", els)
        els <- gsub("Haemosporida", "Haemospororida", els)
        # in this case we can look at unieuk taxonomy
        # P10K oomycotas can be assigned to the eukprot group
        if ("Oomycota" %in% els){
            new_taxonomy <- add_row(new_taxonomy, sp = mnemo,
                                    d = "d__Eukaryota", p = "p__Diaphoretickes",
                                    c = "c__Stramenopiles", o = "o__other_Gyrista",
                                    f = "f__Peronosporomycetes",
                                    g = paste0("g__", genus),
                                    s = paste0("s__", species_name))
            next
        }
        # Heterolobosea sp need some special love
        if (grepl("Heterolobosea", species_name)){
            new_taxonomy <- add_row(new_taxonomy, sp = mnemo,
                                    d = "d__Eukaryota", p = "p__Discoba",
                                    c = "c__Heterolobosea", o = "o__Heterolobosea",
                                    f = "f__Heterolobosea",
                                    g = "g__",
                                    s = paste0("s__", species_name))
            next
        }
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
                                        d = paste0("d__Eukaryota"), p = paste0("p__", new_sk),
                                        c = paste0("c__", new_class), o = paste0("o__", new_order),
                                        f = paste0("f__", new_family), g = paste0("g__", genus),
                                        s = paste0("s__", species_name))
                break
            }
        }
        if (found==FALSE){
            print(species_name)
        }
    }
}

# Arrange the new taxonomy table by taxonomic levels and write it to the output file
new_taxonomy <- arrange(new_taxonomy, p, c, f, o, g) %>% 
    unite(taxon,d,p,c,o,f,g,s, sep = ";")

# write_delim(new_taxonomy, stdout(), delim = "\t", col_names = FALSE)
write.table(new_taxonomy, quote=FALSE, sep = "\t", row.names = FALSE, 
            col.names = FALSE, snakemake@output[[1]])