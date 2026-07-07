# Remove cross-talk (tag-jumps) using the UNCROSS2 score

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Detect and remove per-cell cross-talk (also called tag-jumps) in an OTU
table using Edgar's UNCROSS2 score (Edgar 2018,
[doi:10.1101/400762](https://doi.org/10.1101/400762) ). Unlike
tail-trimming approaches (e.g. `metabaR::tagjumpslayer()`), the score
scales the expected leak by each OTU's **total** abundance: an abundant
genuine OTU is expected to leak more reads into wells where it is absent
than a rare one.

For every OTU (row) an expected per-sample cross-talk count is computed
as `z = f * N / n`, where `N` is the OTU's total abundance across the
whole object and `n` is the **total number of samples** (Edgar's
spreading model: cross-talk leaks into all wells, mostly the zero/low
ones). Each cell then receives the UNCROSS2 score

\$\$s = \frac{2}{1 + \exp(k \cdot x / z)}\$\$

where `x` is the observed count in that cell. The score tends to `1` as
`x -> 0` (looks like cross-talk) and to `0` as `x >> z` (genuine
signal). Cells with `s >= tmin` are treated as cross-talk. The default
`tmin = 0.1` matches Edgar's "weak" (yellow) cut-off; use `0.4` for
"strong" only.

## Usage

``` r
uncross2_pq(
  physeq,
  f = 0.01,
  tmin = 0.1,
  k = 1,
  method = c("cut", "subtract"),
  return_scores = FALSE,
  clean_phyloseq_object = TRUE,
  verbose = TRUE
)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object with integer counts.

- f:

  (numeric or `"auto"`, default `0.01`) Expected cross-talk rate: the
  fraction of an OTU's reads assumed to be misassigned across the plate.
  Set `f = "auto"` to estimate it de novo with
  [`estimate_xtalk_rate_pq()`](https://adrientaudiere.github.io/tidypq/reference/estimate_xtalk_rate_pq.md);
  this errors if the estimate is undetermined.

- tmin:

  (numeric, default `0.1`) Minimum UNCROSS2 score for a cell to be
  flagged as cross-talk. `0.1` = weak, `0.4` = strong (Edgar 2018).

- k:

  (numeric, default `1`) Steepness of the logistic score, applied as a
  slope on the exponent (`exp(k * x / z)`). `k = 1` reproduces UNCROSS2;
  `k > 1` makes the transition sharper, `k < 1` gentler.

- method:

  (character, default `"cut"`) How flagged cross-talk is handled.

  - `"cut"`: set every flagged cell to `0` (UNCROSS2 behaviour).

  - `"subtract"`: subtract the expected count `z` from **all** cells of
    each OTU and floor negatives at `0` (background subtraction; ignores
    `tmin`).

- return_scores:

  (logical, default `FALSE`) If `TRUE`, return a list with the input
  phyloseq object, the cleaned phyloseq object and the score / flag
  matrices instead of just the phyloseq object (see Value).

- clean_phyloseq_object:

  (logical, default `TRUE`) Whether to clean the resulting object with
  [`MiscMetabar::clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.html)
  to drop emptied taxa/samples.

- verbose:

  (logical, default `TRUE`) Whether to print a one-line summary.

## Value

By default a phyloseq object with flagged cross-talk removed. If
`return_scores = TRUE`, a list with:

- `physeq`: the cleaned phyloseq object.

- `scores`: taxa-by-samples matrix of UNCROSS2 scores.

- `tagjump`: logical taxa-by-samples matrix of flagged cells (always
  relative to `tmin`, even when `method = "subtract"`). Matrices are
  oriented taxa-as-rows regardless of the input orientation.

## Relationship to Edgar's usearch binary

This function faithfully implements the *published* UNCROSS2 score (eq.
3), verified cell-for-cell against an independent reference. It differs
from the usearch v11 `otutab_xtalk` binary in two ways that users should
know:

- **Rate.** `f` is fixed (default `0.01`) unless `f = "auto"`, whereas
  the binary *estimates* it de novo. Use
  [`estimate_xtalk_rate_pq()`](https://adrientaudiere.github.io/tidypq/reference/estimate_xtalk_rate_pq.md)
  (or `f = "auto"`) for a de novo estimate.

- **Aggressiveness.** The binary's `-otutabout` filtering is far more
  aggressive than eq. (3): empirically it zeroes essentially every count
  below the per-OTU mean (`N_i/n`), keeping only the few "big" source
  samples, and removed ~40x more cells than eq. (3) at the same rate on
  data. That behaviour is undocumented and beyond the paper, so this
  function does **not** reproduce the binary's filtered table; it is the
  gentler, paper-faithful filter.

UNCROSS2 is defined for **integer count** tables. The usearch binary
rejects non-integer input (`Invalid count`); this function will process
relative abundances but the score's interpretation no longer holds, so
prefer raw counts.

## References

Edgar, R.C. (2018) UNCROSS2: identification of cross-talk in 16S rRNA
OTU tables. *bioRxiv* 400762.
[doi:10.1101/400762](https://doi.org/10.1101/400762)

## See also

[`estimate_xtalk_rate_pq()`](https://adrientaudiere.github.io/tidypq/reference/estimate_xtalk_rate_pq.md)
for the de novo rate;
[`decontam_control_samples_pq()`](https://adrientaudiere.github.io/tidypq/reference/decontam_control_samples_pq.md)
for blank-based decontamination.

## Author

Adrien Taudière
