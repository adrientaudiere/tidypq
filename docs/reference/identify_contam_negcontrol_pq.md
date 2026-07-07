# Detect contaminants from negative-control occurrence patterns

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Examines taxa detected in negative-control samples and assigns each a
contamination sub-type from read abundance and occurrence patterns:

- **artifact**: very low total reads AND present in very few samples –
  likely sequencing noise or cross-talk.

- **lab_contaminant**: predominantly found in negative controls relative
  to real samples – likely introduced during library preparation.

- **sample_contaminant**: detected in controls but also widespread in
  real samples – may represent genuine biological signal.

This is a pure detector: it flags taxa and returns a
[contam_tbl](https://adrientaudiere.github.io/tidypq/reference/contam_tbl.md)
without modifying `physeq`. Because `sample_contaminant` may be genuine
signal, only `artifact` and `lab_contaminant` are flagged by default
(see `flag_categories`). Use
[`filter_contam_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_contam_pq.md)
to remove the flagged taxa.

## Usage

``` r
identify_contam_negcontrol_pq(
  physeq,
  neg_control,
  min_reads_artifact = 10,
  max_samples_artifact = 2,
  max_ratio_lab_contam = 0.2,
  min_neg_samples_lab = 2,
  flag_categories = c("artifact", "lab_contaminant"),
  verbose = TRUE
)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- neg_control:

  An expression evaluated on `sample_data` returning `TRUE` for
  negative-control samples (e.g. `is_control == TRUE`). Use `.` to refer
  to the phyloseq object.

- min_reads_artifact:

  (numeric, default `10`) Maximum total reads for a taxon to be
  considered an artifact.

- max_samples_artifact:

  (integer, default `2`) Maximum total number of samples in which an
  artifact taxon may appear.

- max_ratio_lab_contam:

  (numeric, default `0.2`) Maximum ratio of non-control to control
  sample occurrences for the lab-contaminant category.

- min_neg_samples_lab:

  (integer, default `2`) Minimum number of control samples a
  lab-contaminant taxon must appear in.

- flag_categories:

  (character, default `c("artifact", "lab_contaminant")`) Which
  sub-types to include in the returned `contam_tbl`. Add
  `"sample_contaminant"` to flag those as well.

- verbose:

  (logical, default `TRUE`) Print a one-line detection summary.

## Value

A
[contam_tbl](https://adrientaudiere.github.io/tidypq/reference/contam_tbl.md)
with one row per flagged taxon, `method = "negcontrol"`, and the
evidence columns `subtype` (`artifact`, `lab_contaminant`, or
`sample_contaminant`), `total_reads`, `reads_in_neg`,
`reads_in_samples`, `n_neg_samples`, `n_non_neg_samples`, and
`ratio_non_neg_to_neg`. Empty when no eligible taxon is found.

## See also

[`filter_contam_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_contam_pq.md)
to remove flagged taxa.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
data(data_fungi)
# \donttest{
pq <- mutate_samdata_pq(
  data_fungi,
  is_control = sample_sums(.) < sort(sample_sums(.))[4]
)
flagged <- identify_contam_negcontrol_pq(pq, is_control)
#> ℹ negcontrol: flagged 0 taxa (artifact=0, lab_contaminant=0).
flagged
#> ✔ contam_tbl: no flagged taxa.
#> # A tibble: 0 × 9
#> # ℹ 9 variables: taxon <chr>, method <chr>, subtype <chr>, total_reads <dbl>,
#> #   reads_in_neg <dbl>, reads_in_samples <dbl>, n_neg_samples <dbl>,
#> #   n_non_neg_samples <dbl>, ratio_non_neg_to_neg <dbl>
# }
```
