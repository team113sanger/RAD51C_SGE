## Prepare data from da Silva et al 2019 (Scientific Reports)
## DOI: https://doi.org/10.1038/s41598-019-52078-9
## Single guide, genome-wide CRISPR screen with GeCKO 2.0 library
## Wild type HAP1 vs HAP 1 with DNA Ligase IV (LIG4) knocked-out

# Load libraries
library(readxl)
library(tidyverse)

# Set working directory
top_dir <- '/Users/vo1/cancer/projects/hap1_wt_vs_lig4_daSilva_crispr_screen'
setwd(top_dir)

# Download Supplementary Table 1
# Contains normalised counts (3 x replicates per condition)
input_file <- file.path(top_dir, 'daSilva_Supplementary_Table_1.xlsx')
url <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-019-52078-9/MediaObjects/41598_2019_52078_MOESM2_ESM.xlsx'
download.file(url = url, destfile = input_file, quiet = T)

# Read in data from Excel
# First two columns contain library information (guide sequence and gene name)
# Next six columns contain the normalised replicate counts interleved (WT and L4)
# Column seven contains the normalised plasmid counts (library_representation)

norm_counts <- read_excel(input_file, sheet = "Sheet1", col_names = T, range = cell_cols("A:I"))

# Check total number of guides
num_guides <- norm_counts |> pull(`sgRNA sequence`) |> unique() |> length()
print(paste('Number of guides (norm counts):', num_guides))

# Check total number of genes 
num_genes <- norm_counts |> pull(gene) |> unique() |> length()
print(paste('Number of genes (norm counts):', num_genes))

# Get guides per gene distribution
guides_per_gene <- norm_counts |> count(gene)

# There is something odd here
# GeCKO v2 (https://www.addgene.org/pooled-library/zhang-human-gecko-v2/) 
# Has 19,050 genes across 123,411 guides (~2k controls, 1000 per half library)
# da Silva paper "using the GeCKO v2.0 library that targets 19,052 genes with 122,417 sgRNAs"
# We need to compare the library with the normalised counts to see what's missing

# Download the GeCKO v2 library (guide A)
gecko_A_file <- file.path(top_dir, 'GeCKO_v2_A_gRNA.csv')
gecko_A_url <- 'https://media.addgene.org/cms/filer_public/a4/b8/a4b8d181-c489-4dd7-823a-fe267fd7b277/human_geckov2_library_a_09mar2015.csv'
download.file(url = gecko_A_url, destfile = gecko_A_file, quiet = T)

# Download the GeCKO v2 library (guide B)
gecko_B_file <- file.path(top_dir, 'GeCKO_v2_B_gRNA.csv')
gecko_B_url <- 'https://media.addgene.org/cms/filer_public/2d/8b/2d8baa42-f5c8-4b63-9c6c-bd98f333b29e/human_geckov2_library_b_09mar2015.csv'
download.file(url = gecko_B_url, destfile = gecko_B_file, quiet = T)

# Read in GeCKO v2 libraries and combine
gecko_A <- read_csv(gecko_A_file)
gecko_B <- read_csv(gecko_B_file)
gecko <- rbind(gecko_A[,1:3], gecko_B[,1:3])

# Check number of guides
num_gecko_guides <- gecko |> nrow()
print(paste('Number of guides (GeCKO):', num_gecko_guides))

# Check number of controls
num_gecko_controls <- gecko |> filter(grepl('NonTargetingControlGuideForHuman', gene_id)) |> pull(UID) |> unique() |> length()
print(paste('Number of controls (GeCKO):', num_gecko_controls))

# Check number of genes (minus controls)
num_gecko_genes <- gecko |> filter(!grepl('Control', gene_id)) |> pull(gene_id) |> unique() |> length()
print(paste('Number of genes (GeCKO):', num_gecko_genes))

# Combine library and counts
combined <- full_join(mageck_counts, gecko, by = c('sgRNA' = 'seq')) |> 
  mutate('gene_is_match' = ifelse(gene == gene_id & !is.na(gene) & !is.na(gene_id), 1, 0))

# Fix where the gene name has been parsed as a date in Excel
genes_as_dates <- combined |> filter(grepl('-', gene) & gene_is_match == 0)
combined <- combined |>
  rowwise() |>
  mutate('gene' = ifelse(UID %in% genes_as_dates$UID, gene_id, gene)) |>
  mutate('gene_is_match' = ifelse(gene == gene_id & !is.na(gene) & !is.na(gene_id), 1, 0))

# Get missing guides
missing_guides <- combined |> filter(gene_is_match == 0)
missing_guides_minus_controls <- missing_guides |> filter(!grepl('NonTargetingControlGuideForHuman', gene_id))

# Format counts for use with MAGeCK
mageck_counts <- norm_counts |>
  select('sgRNA' = `sgRNA sequence`, gene, 
         'plasmid' = library_representation, 
         WT_rep1, WT_rep2, WT_rep3, 
         'LIG4_rep1' = L4_rep1, 'LIG4_rep2' = L4_rep2, LIG4_rep3)

# Write counts out to a file
write.table(mageck_counts, 
            file.path(top_dir, 'daSilva_normalised_counts_for_mageck.tsv'),
            row.names = F, quote = F, sep = "\t")