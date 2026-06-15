# Flag taxa whose genus belongs to a known-contaminant blocklist

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Flag taxa whose genus matches a blocklist of known reagent and
laboratory contaminants. Unlike
[`contam_corr_pq()`](https://adrientaudiere.github.io/tidypq/reference/contam_corr_pq.md)
(which infers contaminants from the reads-vs-depth relationship) and
[`neg_control_classify_pq()`](https://adrientaudiere.github.io/tidypq/reference/neg_control_classify_pq.md)
(which uses negative-control occurrence patterns), this method relies
solely on taxonomy and therefore needs no negative controls. It is best
used as a complementary, curation-oriented filter – review flagged taxa
before removal, as several blocklist genera also occur as genuine
environmental taxa.

## Usage

``` r
contam_blocklist_pq(
  physeq,
  genus_blocklist = known_contaminant_genera(),
  tax_rank = "Genus",
  extra_genera = NULL,
  ignore_case = TRUE,
  verbose = TRUE
)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object with a `tax_table`.

- genus_blocklist:

  (character, default:
  [`known_contaminant_genera()`](https://adrientaudiere.github.io/tidypq/reference/known_contaminant_genera.md))
  Vector of genus names to flag.

- tax_rank:

  (character, default: "Genus") Name of the `tax_table` column holding
  genus-level assignments.

- extra_genera:

  (character, default: NULL) Additional genus names to append to
  `genus_blocklist`.

- ignore_case:

  (logical, default: TRUE) If TRUE, matching is case-insensitive.

- verbose:

  (logical, default: TRUE) If TRUE, print a `cli` summary of the flagged
  taxa.

## Value

A named list with four elements:

- `stats`:

  A `data.frame` with one row per taxon: `taxon`, `genus`,
  `is_known_contaminant`, `total_reads`, `prevalence`.

- `contaminants`:

  The subset of `stats` flagged as contaminants.

- `blocklist`:

  The genus vector actually used.

- `summary`:

  A named numeric vector summarizing the result.

## See also

[`known_contaminant_genera()`](https://adrientaudiere.github.io/tidypq/reference/known_contaminant_genera.md),
[`contam_corr_pq()`](https://adrientaudiere.github.io/tidypq/reference/contam_corr_pq.md),
[`neg_control_classify_pq()`](https://adrientaudiere.github.io/tidypq/reference/neg_control_classify_pq.md)

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
data(data_fungi)
# \donttest{
# data_fungi is fungal, so the default bacterial blocklist matches
# nothing; extend it with a genus present in the data for illustration.
res <- contam_blocklist_pq(data_fungi, extra_genera = "Mortierella")
#> ✔ Known-contaminant blocklist summary:
#> • Blocklist genera: 84
#> • Taxa flagged: 1 of 1420 (0.1%)
#> • Reads in flagged taxa: 118
res$summary
#>         taxa_total     n_contaminants   pct_contaminants reads_contaminants 
#>             1420.0                1.0                0.1              118.0 
#>     blocklist_size 
#>               84.0 
head(res$contaminants)
#>     taxon       genus is_known_contaminant total_reads prevalence
#> 1 ASV1215 Mortierella                 TRUE         118          1
# }
```
