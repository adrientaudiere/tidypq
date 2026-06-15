# Detect contaminants by correlating taxon reads with total sample depth

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Identifies potential contaminant taxa by examining the correlation
between each taxon's read count (or relative abundance) and the total
number of reads per sample. The method is inspired by the GRIMER tool
and relies on the observation that contaminant taxa, which are
introduced at a roughly constant level during sample processing, show a
**negative** correlation between their relative abundance and total
sample depth: their proportion is inflated in low-biomass samples.

Conversely, true biological taxa tend to scale with total sample depth,
showing a **positive** correlation between absolute read count and total
reads.

The function computes, for each taxon:

- Pearson correlation between absolute read count and total sample reads
  (`cor_abs`)

- Pearson correlation between relative abundance and total sample reads
  (`cor_rel`)

- Slope and p-value from a linear model of relative abundance ~ total
  reads

Taxa are flagged as potential contaminants when `cor_rel` is below a
configurable threshold (default: -0.5).

## Usage

``` r
contam_corr_pq(
  physeq,
  min_total_reads = 1000,
  min_prevalence = 2,
  contam_threshold = -0.5,
  tax_rank = "Species",
  top_n = 20,
  plot = TRUE,
  verbose = TRUE
)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- min_total_reads:

  (numeric, default `1000`) Samples with fewer total reads are excluded
  from the analysis to avoid spurious correlations driven by extreme
  low-depth outliers.

- min_prevalence:

  (integer, default `2`) Taxa present in fewer than this many samples
  are excluded from the analysis.

- contam_threshold:

  (numeric, default `-0.5`) Taxa with `cor_rel` below this value are
  flagged as potential contaminants. Use values closer to 0 for a more
  conservative (fewer flagged) approach, or more negative for a more
  aggressive detection.

- tax_rank:

  (character, default `"Species"`) The taxonomic rank column used for
  labelling taxa in the plot. If the column does not exist, taxa names
  are used instead.

- top_n:

  (integer, default `20`) Number of top contaminant and top
  non-contaminant taxa to highlight in the diagnostic plot.

- plot:

  (logical, default `TRUE`) Whether to produce the diagnostic plot.

- verbose:

  (logical, default `TRUE`) Print a summary of detected contaminants.

## Value

A list with the following elements:

- `stats`: A data.frame with one row per taxon, columns: `taxon`,
  `total_reads`, `prevalence`, `cor_abs`, `cor_rel`, `slope`, `p_value`,
  `is_contaminant`.

- `contaminants`: Subset of `stats` for taxa flagged as contaminants.

- `plot`: A ggplot object (if `plot = TRUE`) with three panels: (A)
  histogram of `cor_rel` values, (B) scatter plot of `cor_abs` vs
  `cor_rel`, (C) read-count vs total-depth for the top flagged taxa.

## See also

[`neg_control_diag_pq()`](https://adrientaudiere.github.io/tidypq/reference/neg_control_diag_pq.md),
[`decontam_sam_control()`](https://adrientaudiere.github.io/tidypq/reference/decontam_sam_control.md)

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
data(data_fungi)

# Inject a fake contaminant (constant reads across samples)
# so the example produces visible contaminant detection results
otu <- as(phyloseq::otu_table(data_fungi), "matrix")
fake <- matrix(rep(50, nrow(otu)), ncol = 1,
  dimnames = list(rownames(otu), "fake_contaminant"))
otu <- cbind(otu, fake)
tax <- rbind(
  as.matrix(phyloseq::tax_table(data_fungi)),
  fake_contaminant = rep("Contaminant", ncol(phyloseq::tax_table(data_fungi)))
)
data_fungi_c <- phyloseq::phyloseq(
  phyloseq::otu_table(otu, taxa_are_rows = FALSE),
  phyloseq::tax_table(tax),
  phyloseq::sample_data(data_fungi)
)

# Detect contaminants
res <- contam_corr_pq(data_fungi_c)
#> Warning: the standard deviation is zero
#> ✔ Contaminant detection summary:
#> • Samples analyzed: 143 (of 185)
#> • Taxa analyzed: 1200
#> • Potential contaminants (cor_rel < -0.5): 1 (0.1%)
#> • Median cor_rel for contaminants: -0.667
#> • Median cor_rel for non-contaminants: -0.001
#> Error: Can't find method for generic `|(e1, e2)`:
#> - e1: <ggplot2::ggplot>
#> - e2: <ggplot2::ggplot>
res$summary
#> Error: object 'res' not found
res$plot
#> Error: object 'res' not found
head(res$contaminants)
#> Error: object 'res' not found

if (FALSE) { # \dontrun{
# Adjust threshold
 res_strict <- contam_corr_pq(data_fungi_c, contam_threshold = -0.2, tax_rank = "Genus_species")
 res$plot
} # }
```
