# Plot UNCROSS2 cross-talk scores against OTU abundance

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Diagnostic scatter plot for the output of
[`uncross2_pq()`](https://adrientaudiere.github.io/tidypq/reference/uncross2_pq.md)
run with `return_scores = TRUE`. Each point is one non-zero occurrence
(an OTU in a sample), placed by the OTU's total abundance across the
object (x) and the read count of that occurrence (y), and coloured by
its UNCROSS2 score. Occurrences scoring at or below `tmin` (kept as
genuine signal) are drawn in grey; those at or above `tmin` (flagged as
cross-talk) are coloured on the score scale. The grey 1:1 line marks
where an occurrence equals the OTU total: genuine "source" samples sit
near it, cross-talk sits far below. Both axes are on a log10 scale.

## Usage

``` r
plot_uncross2_scores(scores_list, grey_non_significant = TRUE)
```

## Arguments

- scores_list:

  (list, required) The list returned by
  [`uncross2_pq()`](https://adrientaudiere.github.io/tidypq/reference/uncross2_pq.md)
  with `return_scores = TRUE`. Must contain `old_physeq` (the input
  object), the `scores` matrix, and the `tmin` threshold.

- grey_non_significant:

  (logical, default `TRUE`) Draw occurrences scoring at or below `tmin`
  in grey rather than on the colour scale.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## See also

[`uncross2_pq()`](https://adrientaudiere.github.io/tidypq/reference/uncross2_pq.md),
which produces the scored list.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
data(data_fungi)
# \donttest{
res <- uncross2_pq(data_fungi, return_scores = TRUE)
#> ℹ uncross2 (cut): removed 984 reads across 250620 flagged cells (occurences) (f = 0.01, tmin = 0.1, k = 1).
plot_uncross2_scores(res)
#> Taxa are now in rows.

# }
```
