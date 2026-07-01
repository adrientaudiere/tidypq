# Convert a phyloseq object to a tidy tibble

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

`pq_to_tidy()` is the canonical phyloseq-to-tidy data-preparation module
for the pqverse. It owns the full pipeline: cleaning, flattening,
optional sample aggregation, optional abundance transformation, and
fact/bifactor resolution. Every ggplotpq plot function calls
`pq_to_tidy()` instead of reinventing its own data-prep pipeline.

The function builds the tidy tibble from scratch using
[`tidyr::pivot_longer()`](https://tidyr.tidyverse.org/reference/pivot_longer.html)
(it does not wrap
[`phyloseq::psmelt()`](https://rdrr.io/pkg/phyloseq/man/psmelt.html)).
Taxonomic rank grouping happens in tibble space (no
[`tax_glom()`](https://rdrr.io/pkg/phyloseq/man/tax_glom.html)). Sample
aggregation happens in tibble space (no
[`merge_samples()`](https://rdrr.io/pkg/phyloseq/man/merge_samples-methods.html)
round-trip).

**Alternative transform approaches:** users may also transform the
phyloseq object before calling `pq_to_tidy()` via
[`MiscMetabar::transform_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/transform_pq.html),
or transform the tibble after via tidyverse
[`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html) (e.g.
`dplyr::mutate(abundance = abundance / sum(abundance))`).

## Usage

``` r
pq_to_tidy(
  physeq,
  fact = NULL,
  bifactor = NULL,
  merge_sample_by = NULL,
  transform = NULL,
  ranks = NULL,
  filter_zero = TRUE,
  verbose = TRUE
)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object. Must have `sample_data`.

- fact:

  (character or NULL, default NULL) Name of a `sample_data` column used
  as the primary grouping factor. `NULL` defaults to `sample_id`.

- bifactor:

  (character or NULL, default NULL) Name of a `sample_data` column used
  as the secondary grouping factor. Must have exactly 2 levels. `NULL`
  defaults to `sample_id`.

- merge_sample_by:

  (character or NULL, default NULL) Name of a `sample_data` column used
  to aggregate (merge) samples. Abundances are summed per group; other
  columns take their first value. `sample_id` becomes the group label.

- transform:

  (string, function, or NULL, default NULL) Abundance transformation
  applied per `sample_id` after aggregation. A string is resolved via
  [`match.fun()`](https://rdrr.io/r/base/match.fun.html) (e.g.
  `"log1p"`, `"sqrt"`). A function is applied directly. `NULL` leaves
  raw counts.

- ranks:

  (character or NULL, default NULL) Taxonomic ranks to include as
  columns. `NULL` includes all ranks from `tax_table`. A character
  vector (e.g. `c("Phylum", "Genus")`) selects specific ranks. Ignored
  if no `tax_table`.

- filter_zero:

  (logical, default TRUE) Drop rows where `abundance == 0`.

- verbose:

  (logical, default TRUE) Print diagnostic messages.

## Value

An ungrouped tibble with columns:

- `sample_id` (character) — sample identifier or merge group label.

- `taxon_id` (character) — OTU/ASV identifier.

- `abundance` (numeric) — counts or transformed values.

- `abundance_raw` (numeric) — pre-transform counts, present only when
  `transform` is applied.

- `fact` (character) — resolved primary grouping column.

- `bifactor` (character) — resolved secondary grouping column.

- Rank columns (character) — one per requested rank, `NA` replaced by
  `"Unknown"`. Present only if `tax_table` exists.

- All other `sample_data` columns (original names and types).

## See also

[`phyloseq::psmelt()`](https://rdrr.io/pkg/phyloseq/man/psmelt.html),
[`MiscMetabar::psmelt_samples_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/psmelt_samples_pq.html),
[`MiscMetabar::clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.html),
[`MiscMetabar::transform_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/transform_pq.html)

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
# Basic tidy conversion
tidy_df <- pq_to_tidy(data_fungi_mini)
#> After clean_pq: 137 samples, 45 taxa.
#> Filtered 5605 zero-abundance rows.
#> Returning tibble with 560 rows and 24 columns.
head(tidy_df)
#> # A tibble: 6 × 24
#>   taxon_id sample_id     abundance X     Sample_names Tree_name Sample_id Height
#>   <chr>    <chr>             <dbl> <chr> <chr>        <chr>         <int> <chr> 
#> 1 ASV7     A10-005-B_S1…       767 A10-… A10-005-B_S… A10-005         188 Low   
#> 2 ASV7     A10-005-H_S1…     10187 A10-… A10-005-H_S… A10-005         189 High  
#> 3 ASV7     A10-005-M_S1…       263 A10-… A10-005-M_S… A10-005         190 Middle
#> 4 ASV7     A12-007_S191…        35 A12-… A12-007_S191 A12-007         191 NA    
#> 5 ASV7     A8-005_S4_ME…      4215 A8-0… A8-005_S4    A8-005            4 NA    
#> 6 ASV7     AD26-005-M_S…         1 AD26… AD26-005-M_… AD26-005         11 Middle
#> # ℹ 16 more variables: Diameter <chr>, Time <int>, Domain <chr>, Phylum <chr>,
#> #   Class <chr>, Order <chr>, Family <chr>, Genus <chr>, Species <chr>,
#> #   Trophic.Mode <chr>, Guild <chr>, Trait <chr>, Confidence.Ranking <chr>,
#> #   Genus_species <chr>, fact <chr>, bifactor <chr>

# With fact and bifactor
tidy_df <- pq_to_tidy(data_fungi_mini, fact = "Height", bifactor = "Time")
#> After clean_pq: 137 samples, 45 taxa.
#> Error in pq_to_tidy(data_fungi_mini, fact = "Height", bifactor = "Time"): bifactor 'Time' must have exactly 2 levels, found 5.

# With sample aggregation
tidy_df <- pq_to_tidy(data_fungi_mini, merge_sample_by = "Height")
#> After clean_pq: 137 samples, 45 taxa.
#> Aggregating samples by 'Height'.
#> Filtered 38 zero-abundance rows.
#> Returning tibble with 142 rows and 23 columns.

# With log1p transform
tidy_df <- pq_to_tidy(data_fungi_mini, transform = "log1p")
#> After clean_pq: 137 samples, 45 taxa.
#> Applying transform.
#> Filtered 5605 zero-abundance rows.
#> Returning tibble with 560 rows and 25 columns.

# With selected ranks only
tidy_df <- pq_to_tidy(data_fungi_mini, ranks = c("Phylum", "Genus"))
#> After clean_pq: 137 samples, 45 taxa.
#> Filtered 5605 zero-abundance rows.
#> Returning tibble with 560 rows and 14 columns.
```
