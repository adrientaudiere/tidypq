# Detect contaminants by correlating taxon reads with total sample depth

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Flags potential contaminant taxa from the correlation between each
taxon's relative abundance and total sample depth. Contaminants,
introduced at a roughly constant level during processing, show a
**negative** correlation (their proportion is inflated in low-biomass
samples). Inspired by the GRIMER tool.

This is a pure detector: it flags taxa and returns a
[contam_tbl](https://adrientaudiere.github.io/tidypq/reference/contam_tbl.md)
without modifying `physeq`. Use
[`filter_contam_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_contam_pq.md)
to remove the flagged taxa and
[`plot.contam_tbl()`](https://adrientaudiere.github.io/tidypq/reference/autoplot.contam_tbl.md)
for a quick diagnostic.

## Usage

``` r
identify_contam_corr_pq(
  physeq,
  min_total_reads = 1000,
  min_prevalence = 2,
  contam_threshold = -0.5,
  verbose = TRUE
)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- min_total_reads:

  (numeric, default `1000`) Samples with fewer total reads are excluded
  to avoid spurious correlations from low-depth outliers.

- min_prevalence:

  (integer, default `2`) Taxa present in fewer than this many samples
  are excluded from the analysis.

- contam_threshold:

  (numeric, default `-0.5`) Taxa whose relative- abundance correlation
  (`cor_rel`) is below this value are flagged. Use values closer to 0 to
  flag more aggressively.

- verbose:

  (logical, default `TRUE`) Print a one-line detection summary.

## Value

A
[contam_tbl](https://adrientaudiere.github.io/tidypq/reference/contam_tbl.md)
with one row per flagged taxon, `method = "corr"`, and the evidence
columns `total_reads`, `prevalence`, `cor_abs`, `cor_rel`, `slope`, and
`p_value`. Empty when nothing is flagged.

## See also

[`filter_contam_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_contam_pq.md)
to remove flagged taxa,
[contam_tbl](https://adrientaudiere.github.io/tidypq/reference/contam_tbl.md)
for the return type.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
# \donttest{
# Inject a constant-reads fake contaminant so the example flags something
otu <- as(phyloseq::otu_table(data_fungi), "matrix")
fake <- matrix(rep(50, nrow(otu)), ncol = 1,
  dimnames = list(rownames(otu), "fake_contaminant"))
tax <- rbind(
  as.matrix(phyloseq::tax_table(data_fungi)),
  fake_contaminant = rep("Contaminant", ncol(phyloseq::tax_table(data_fungi)))
)
data_fungi_c <- phyloseq::phyloseq(
  phyloseq::otu_table(cbind(otu, fake), taxa_are_rows = FALSE),
  phyloseq::tax_table(tax),
  phyloseq::sample_data(data_fungi)
)

flagged <- identify_contam_corr_pq(data_fungi_c)
#> ℹ corr: flagged 1 of 1200 taxa (0.1%) at cor_rel < -0.5.
flagged
#> contam_tbl: 1 flagged taxon from method "corr"
#> # A tibble: 1 × 8
#>   taxon          method total_reads prevalence cor_abs cor_rel    slope  p_value
#>   <chr>          <chr>        <dbl>      <dbl>   <dbl>   <dbl>    <dbl>    <dbl>
#> 1 fake_contamin… corr          7150        143      NA  -0.667 -6.32e-7 8.96e-20
# }
```
