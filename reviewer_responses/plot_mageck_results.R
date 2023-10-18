## Compare MAGeCK results generated from normalised counts from da Silva et al 2019 (Scientific Reports)
## DOI: https://doi.org/10.1038/s41598-019-52078-9
## Single guide, genome-wide CRISPR screen with GeCKO 2.0 library
## Wild type HAP1 vs HAP 1 with DNA Ligase IV (LIG4) knocked-out

# Load libraries
library(tidyverse)
library(scales)
library(ggpubr)
library(ggsci)
library(ggrepel)

# Set working directory
top_dir <- '/Users/vo1/cancer/projects/hap1_wt_vs_lig4_daSilva_crispr_screen'
setwd(top_dir)

# Read in results for MAGeCK contrasts
wt_vs_plasmid <- read.table(file.path(top_dir, 'MAGeCK', 'WT_vs_plasmid.gene_summary.txt'), header = T, sep = "\t")
lig4_vs_plasmid <- read.table(file.path(top_dir, 'MAGeCK', 'LIG4_vs_plasmid.gene_summary.txt'), header = T, sep = "\t")
lig4_vs_wt <- read.table(file.path(top_dir, 'MAGeCK', 'LIG4_vs_WT.gene_summary.txt'), header = T, sep = "\t")

# Label contrasts
lig4_vs_plasmid <- lig4_vs_plasmid |> rename_at(vars(-id), function( x ) paste('LIG4vsPlasmid' , x, sep = '.'))
wt_vs_plasmid <- wt_vs_plasmid |> rename_at(vars(-id), function( x ) paste('WTvsPlasmid' , x, sep = '.'))
lig4_vs_wt <- lig4_vs_wt |> rename_at(vars(-id), function( x ) paste('LIG4vsWT' , x, sep = '.'))

# Merge contrasts
combined_data <- wt_vs_plasmid |>
  full_join(lig4_vs_plasmid, by = 'id') |> 
  full_join(lig4_vs_wt, by = 'id')
  
# Classify results
fdr_cutoff <- 0.05
combined_data <- combined_data |>
  mutate('WTvsPlasmid.enriched' = ifelse(WTvsPlasmid.pos.fdr <= fdr_cutoff & WTvsPlasmid.pos.lfc > 0, TRUE, FALSE)) |> 
  mutate('WTvsPlasmid.depleted' = ifelse(WTvsPlasmid.neg.fdr <= fdr_cutoff & WTvsPlasmid.neg.lfc < 0, TRUE, FALSE)) |>
  mutate('LIG4vsPlasmid.enriched' = ifelse(LIG4vsPlasmid.pos.fdr <= fdr_cutoff & LIG4vsPlasmid.pos.lfc > 0, TRUE, FALSE)) |> 
  mutate('LIG4vsPlasmid.depleted' = ifelse(LIG4vsPlasmid.neg.fdr <= fdr_cutoff & LIG4vsPlasmid.neg.lfc < 0, TRUE, FALSE)) |>
  mutate('LIG4vsWT.enriched' = ifelse(LIG4vsWT.pos.fdr <= fdr_cutoff & LIG4vsWT.pos.lfc > 0, TRUE, FALSE)) |> 
  mutate('LIG4vsWT.depleted' = ifelse(LIG4vsWT.neg.fdr <= fdr_cutoff & LIG4vsWT.neg.lfc < 0, TRUE, FALSE)) |> 
  mutate('sig.enriched' = ifelse(LIG4vsWT.enriched == TRUE & LIG4vsPlasmid.enriched == TRUE, TRUE, FALSE)) |>
  mutate('sig.depleted'  = ifelse(LIG4vsWT.depleted == TRUE & LIG4vsPlasmid.depleted == TRUE, TRUE, FALSE)) |>
  mutate('sig'  =  ifelse(sig.enriched == TRUE | sig.depleted == TRUE , TRUE, FALSE))

# Save results
write.table(combined_data, file.path(top_dir, 'mageck_results.tsv'), sep = "\t", row.names = F, quote = F)

# Create plot data
plot_data <- combined_data |>
  mutate(direction = ifelse((LIG4vsWT.enriched & !LIG4vsPlasmid.enriched) | (LIG4vsWT.enriched & LIG4vsPlasmid.enriched), 'enriched',
                     ifelse((LIG4vsWT.depleted & !LIG4vsPlasmid.depleted) | (LIG4vsWT.depleted & LIG4vsPlasmid.depleted), 'depleted', NA)))

# Specify the levels for legends
plot_data$direction <- factor(plot_data$direction, levels = c("depleted", "enriched"))
plot_data$sig <- factor(plot_data$sig, levels = c("TRUE","FALSE"))

# Plot results
genes_to_highlight <- c('BRCA1', 'BRCA2', 'ATM', 
                        'BARD1', 'BRIP1', 'CHEK2', 'NBS1',
                        'NBN', 'PALB2', 'RAD51C', 'RAD51D')
mageck_plot <-
  ggplot(plot_data, aes(x = LIG4vsPlasmid.neg.lfc, y = LIG4vsWT.neg.lfc, label = id, shape = sig, color = direction)) +
    geom_point(data = subset(plot_data, is.na(direction) == T), color = 'gray', alpha = 0.3) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point(data = subset(plot_data, id %in% genes_to_highlight), color = 'forestgreen', alpha = 0.7, size = 3) +
    geom_label_repel(data = subset(plot_data, id %in% genes_to_highlight), color = 'forestgreen', force = 10, nudge_x = 0.2, max.overlaps = 20) +
    geom_point(data = subset(plot_data, (direction == 'enriched' | direction == 'depleted') & sig == T), alpha = 1, size = 3) +
    geom_label_repel(data = subset(plot_data, direction == 'depleted' & sig == T), color = 'firebrick', force = 10, nudge_x = -0.2, max.overlaps = 20) +
    scale_x_continuous(breaks = pretty_breaks(10), limits = c(-3.5, 1.5), name = 'LIG4 KO vs Plasmid') +
    scale_y_continuous(breaks = pretty_breaks(8), limits = c(-2.5, 2), name = 'LIG4 KO vs WT') +
    scale_shape_discrete(name = 'Significant (FDR < 0.05)') +
    scale_color_manual(name = 'Direction', values = c('enriched' = 'navyblue', 'depleted' = 'firebrick')) +
    theme_pubr(base_size = 18)

# Save plot
ggsave(file.path(top_dir, 'mageck_results.png'), mageck_plot, device = 'png', dpi = 200, width = 10, height = 10)
