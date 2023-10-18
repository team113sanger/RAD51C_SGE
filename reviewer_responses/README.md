# HAP1 WT vs LIG4 KO single guide CRISPR screen

Normalised counts downloaded for a single guide, genome-wide CRISPR screen comparing wild type HAP1 cells with a LIG4 knock out from:

> Ferreira da Silva, J., Salic, S., Wiedner, M. et al. 
> Genome-scale CRISPR screens are efficient in non-homologous end-joining deficient cells. 
> Sci Rep 9, 15751 (2019). 
> [https://doi.org/10.1038/s41598-019-52078-9](https://doi.org/10.1038/s41598-019-52078-9)

Genome-wide screen performed with the human GeCKO v2 library:

> Sanjana, N., Shalem, O. & Zhang, F. 
> Improved vectors and genome-wide libraries for CRISPR screening. 
> Nat Methods 11, 783–784 (2014). 
> [https://doi.org/10.1038/nmeth.3047](https://doi.org/10.1038/nmeth.3047)

Library can be downloaded from:

[https://www.addgene.org/pooled-library/zhang-human-gecko-v2/](https://www.addgene.org/pooled-library/zhang-human-gecko-v2/)

## Download and prepare normalised counts 

`download_daSilva_norm_counts.R ` was used to download and prepare data for analysis with [MAGeCK v0.5.9.3](https://sourceforge.net/p/mageck/wiki/Home/).

* 3 x wild type HAP1 samples 
* 3 x HAP1 LIG4 knock out samples
* 1 x plasmid

Note: some preprocessing has been performed as normalised counts are not present for all GeCKO v2 guides

| File(s) | Description |
| --- | --- |
| `GeCKO_v2_A_gRNA.csv ` and `GeCKO_v2_B_gRNA.csv` | GeCKO v2 CRISPR library from AddGene |
| `daSilva_Supplementary_Table_1.xlsx ` | Supplementary Table 1 with normalised counts |
| `SuppTable1_excl_guides.tsv` | Guides missing from normalised counts |
| `daSilva_normalised_counts_for_mageck.tsv` | Normalised counts for use with MAGeCK |

## MAGeCK

MAGeCK v0.5.9.3](https://sourceforge.net/p/mageck/wiki/Home/) was used to determine enriched and depleted genes for three contrasts using the count matrix prepared in the previous step (`daSilva_normalised_counts_for_mageck.tsv`):

* `WT` (3 replicates) vs `Plasmid` 
* `LIG4` (3 replicates) vs `Plasmid`
* `LIG4` (3 replicates) vs `WT` (3 replicates)

Jobs were run using on a HPC using LSF (`mageck_lsf_jobscript.sh`). Raw results can be found in the `MAGeCK` subdirectory. 

## Post-analysis

MAGeCK results were combined and plotted using `plot_mageck_results.R`. 

| Column | Condition | Description |
| --- | --- | --- |
| `[contrast].enriched` | `TRUE` if `[contrast].pos.fdr` < 0.05 and `[contrast].pos.lfc` > 0 | Whether a gene is significantly enriched in that particular contrast |
| `[contrast].depleted` | `TRUE` if `[contrast].neg.fdr` < 0.05 and `[contrast].neg.lfc` < 0 | Whether a gene is significantly depleted in that particular contrast |
| `sig.enriched` | `TRUE` if `LIG4vsWT.enriched` == `TRUE` & `LIG4vsPlasmid.enriched` == `TRUE` | Whether a gene is enriched in LIG4 KO when compared to both wild type and plasmid |
| `sig.depleted` | `TRUE` if | Whether a gene is depleted in LIG4 KO when compared to both wild type and plasmid |
| `sig` | `TRUE` if `sig.enriched` == `TRUE` or `sig.depleted` == `TRUE` | Whether gene is either significantly enriched or depleted in LIG4 KO when compared to both wild type and plasmid |

Results are summarised in `mageck_results.png` and `mageck_results.tsv`.