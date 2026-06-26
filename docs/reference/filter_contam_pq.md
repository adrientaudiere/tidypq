# Remove contaminant taxa flagged by a detector

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

The single removal verb for contamination handling. It consumes any
`contam_tbl` produced by an `identify_contam_*_pq()` detector (or
several, combined with [`rbind()`](https://rdrr.io/r/base/cbind.html))
and prunes the flagged taxa from `physeq` by name. When the `contam_tbl`
is empty, `physeq` is returned unchanged.

## Usage

``` r
filter_contam_pq(physeq, contam, clean_phyloseq_object = TRUE, verbose = TRUE)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- contam:

  (contam_tbl, required) Output of an `identify_contam_*_pq()` detector.
  Removal is keyed on the `taxon` column, so taxa are matched by name
  rather than position.

- clean_phyloseq_object:

  (logical, default `TRUE`) Whether to clean the resulting object with
  [`MiscMetabar::clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.html)
  to drop emptied taxa/samples.

- verbose:

  (logical, default `TRUE`) Whether to report how many taxa were
  removed.

## Value

A phyloseq object with the flagged taxa removed.

## See also

[`identify_contam_corr_pq()`](https://adrientaudiere.github.io/tidypq/reference/identify_contam_corr_pq.md)
and the other `identify_contam_*_pq()` detectors that produce a
`contam_tbl`.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
# \donttest{
flagged <- identify_contam_corr_pq(data_fungi, verbose = FALSE)
data_clean <- filter_contam_pq(data_fungi, flagged)
#> ✔ No contaminant taxa to remove.
# }
```
