# Estimate the cross-talk rate de novo (Edgar 2018)

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Estimate the cross-talk rate `f` used by
[`uncross2_pq()`](https://adrientaudiere.github.io/tidypq/reference/uncross2_pq.md)
directly from the OTU table, without control/mock samples, following the
de novo procedure of Edgar (2018,
[doi:10.1101/400762](https://doi.org/10.1101/400762) ).

For each OTU `i` with total abundance `N_i`, a count `c` is deemed *low*
(i.e. consistent with cross-talk) when `0 < c <= s*N_i/n`, where `n` is
the number of samples (default `s = 0.1`, so "ten times below the
mean"). Let `S_i` be the number of low counts and `X_i` their sum. An
OTU is a **candidate** when `N_i >= nmin` and `S_i >= smin`, and its
rate is

- `rate = "edgar2018"` (paper eq. 2): `f_i = (X_i/N_i)*(n/S_i)`

- `rate = "fraction"`: `f_i = X_i/N_i`, the raw low-read fraction, which
  is what the usearch `otutab_xtalk` report prints as each OTU's rate.

Candidates with `f_i > fmax` are rejected. If fewer than `cmin`
candidates remain the estimate is **undetermined** and `NA` is returned
(a real UNCROSS2 outcome; several datasets in Edgar's Table 5 are
undetermined). Otherwise the returned rate is the **median** `f_i` over
candidates.

## Usage

``` r
estimate_xtalk_rate_pq(
  physeq,
  s = 0.1,
  nmin = 1000,
  smin = 3L,
  fmax = 0.02,
  cmin = 10L,
  rate = c("edgar2018", "fraction"),
  verbose = TRUE
)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object with integer counts.

- s:

  (numeric, default `0.1`) Low-count fraction of the OTU mean.

- nmin:

  (numeric, default `1000`) Minimum OTU total to be a candidate.

- smin:

  (integer, default `3`) Minimum number of low counts to be a candidate.

- fmax:

  (numeric, default `0.02`) Candidates with an estimated rate above this
  plausible maximum are rejected.

- cmin:

  (integer, default `10`) Minimum number of candidates required; below
  this the estimate is undetermined (`NA`).

- rate:

  (character, default `"edgar2018"`) Per-OTU rate formula, see
  Description.

- verbose:

  (logical, default `TRUE`) Print a one-line summary.

## Value

A single numeric: the median de novo cross-talk rate, or `NA_real_` when
undetermined. The number of candidate OTUs is attached as the
`"n_candidates"` attribute.

## Agreement with usearch

This reimplementation follows the *paper*. The usearch v11 binary
diverges from its own description in undocumented ways (an absolute
low-count cap of ~8 rather than `s*N_i/n`, the `X_i/N_i` rate, and an
iterative candidate selection), so the returned value is of the **same
order of magnitude** as usearch's de novo estimate but is **not**
bit-exact. On integer-count data the per-OTU `X_i/N_i` reproduces
usearch's reported rate when the same low-count set is used; the median
over the whole table does not. See
[`uncross2_pq()`](https://adrientaudiere.github.io/tidypq/reference/uncross2_pq.md)
for the score-side comparison.

## References

Edgar, R.C. (2018) UNCROSS2: identification of cross-talk in 16S rRNA
OTU tables. *bioRxiv* 400762.
[doi:10.1101/400762](https://doi.org/10.1101/400762)

## See also

[`uncross2_pq()`](https://adrientaudiere.github.io/tidypq/reference/uncross2_pq.md),
which accepts `f = "auto"` to use this estimate.

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
  estimate_xtalk_rate_pq(data_fungi)
} # }
```
