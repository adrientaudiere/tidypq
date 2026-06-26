# Detect taxa whose genus belongs to a known-contaminant blocklist

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Flags taxa whose genus matches a blocklist of known reagent and
laboratory contaminants. A taxonomy-only method requiring no negative
controls, best used as a complementary, curation-oriented filter –
review flagged taxa before removal, as several blocklist genera also
occur as genuine environmental taxa.

This is a pure detector: it flags taxa and returns a
[contam_tbl](https://adrientaudiere.github.io/tidypq/reference/contam_tbl.md)
without modifying `physeq`. Use
[`filter_contam_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_contam_pq.md)
to remove the flagged taxa.

## Usage

``` r
identify_contam_blocklist_pq(
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

  (character, default
  [`known_contaminant_genera()`](https://adrientaudiere.github.io/tidypq/reference/known_contaminant_genera.md))
  Vector of genus names to flag.

- tax_rank:

  (character, default `"Genus"`) Name of the `tax_table` column holding
  genus-level assignments.

- extra_genera:

  (character, default `NULL`) Additional genus names to append to
  `genus_blocklist`.

- ignore_case:

  (logical, default `TRUE`) If `TRUE`, matching is case-insensitive.

- verbose:

  (logical, default `TRUE`) Print a one-line detection summary.

## Value

A
[contam_tbl](https://adrientaudiere.github.io/tidypq/reference/contam_tbl.md)
with one row per flagged taxon, `method = "blocklist"`, and the evidence
columns `genus`, `total_reads`, and `prevalence`. Empty when no genus
matches the blocklist.

## See also

[`known_contaminant_genera()`](https://adrientaudiere.github.io/tidypq/reference/known_contaminant_genera.md)
for the default blocklist,
[`filter_contam_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_contam_pq.md)
to remove flagged taxa.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
data(data_fungi)
# \donttest{
# data_fungi is fungal, so the default bacterial blocklist matches nothing;
# extend it with a genus present in the data for illustration.
identify_contam_blocklist_pq(data_fungi, extra_genera = "Mortierella")
#> ℹ blocklist: flagged 1 of 1420 taxa from 84 blocklist genera.
#> contam_tbl: 1 flagged taxon from method "blocklist"
#> # A tibble: 1 × 5
#>   taxon   method    genus       total_reads prevalence
#>   <chr>   <chr>     <chr>             <dbl>      <dbl>
#> 1 ASV1215 blocklist Mortierella         118          1
# }
```
