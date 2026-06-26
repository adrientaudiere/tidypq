# The `contam_tbl` class

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

A `contam_tbl` is the canonical object returned by every tidypq
contaminant detector (`identify_contam_*_pq()`). It is a tibble subclass
with one row per flagged taxon and at least two columns:

- `taxon`:

  Taxon name, matching
  [`phyloseq::taxa_names()`](https://rdrr.io/pkg/phyloseq/man/taxa_names-methods.html).

- `method`:

  Short token naming the detector that flagged the taxon (e.g. `"corr"`,
  `"blocklist"`, `"primer"`).

Detectors add their own method-specific evidence columns (correlations,
p-values, matched primers, …). A detector that flags nothing returns a
zero-row `contam_tbl`, never `NULL`, so outputs always compose.

Outputs of several detectors can be row-bound with
[`rbind()`](https://rdrr.io/r/base/cbind.html); evidence columns absent
from one detector are filled with `NA`. The combined table is consumed
by
[`filter_contam_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_contam_pq.md)
(removal) and by
[`plot.contam_tbl()`](https://adrientaudiere.github.io/tidypq/reference/autoplot.contam_tbl.md)
(diagnostics).

## Usage

``` r
new_contam_tbl(x = tibble::tibble(taxon = character(), method = character()))

is_contam_tbl(x)
```

## Arguments

- x:

  (data.frame, required) A data frame carrying at least `taxon` and
  `method` columns.

## Value

`new_contam_tbl()` returns a `contam_tbl`. `is_contam_tbl()` returns a
length-one logical.

## Author

Adrien Taudière
