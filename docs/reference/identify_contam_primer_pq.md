# Detect taxa whose representative sequence contains a primer

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

For each taxon, checks whether the representative sequence in the
`@refseq` slot matches any primer in `primers` (forward or reverse
complement) using IUPAC-aware pattern matching. Taxa embedding a primer
are likely chimeric or otherwise artefactual.

This is a pure detector: it flags taxa and returns a
[contam_tbl](https://adrientaudiere.github.io/tidypq/reference/contam_tbl.md)
without modifying `physeq`. Use
[`filter_contam_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_contam_pq.md)
to remove the flagged taxa.

## Usage

``` r
identify_contam_primer_pq(physeq, primers, verbose = TRUE)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object with a populated `@refseq` slot
  containing DNA sequences.

- primers:

  (named character vector, required) Primer sequences to search for.
  IUPAC ambiguity codes are supported (e.g. `M`, `R`, `Y`, `N`). Names
  are reported in the `matched_primers` evidence column.

- verbose:

  (logical, default `TRUE`) Print a one-line detection summary.

## Value

A
[contam_tbl](https://adrientaudiere.github.io/tidypq/reference/contam_tbl.md)
with one row per flagged taxon, `method = "primer"`, and the evidence
columns `matched_primers` (comma-separated primer names) and `n_reads`.
Empty when no taxon embeds a primer.

## See also

[`filter_contam_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_contam_pq.md)
to remove flagged taxa,
[`MiscMetabar::cutadapt_remove_primers()`](https://adrientaudiere.github.io/MiscMetabar/reference/cutadapt_remove_primers.html)
for cutadapt-based primer removal of fastq files.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
# Requires a phyloseq object with a populated @refseq slot.
primers <- c(
  fwd = "CCCTACGGGGTGCASCAG",
  rev = "GGACTACVSGGGTATCTAAT"
)
identify_contam_primer_pq(data_fungi, primers)
#> ℹ primer: flagged 0 of 1420 taxa embedding a primer.
#> ✔ contam_tbl: no flagged taxa.
#> # A tibble: 0 × 4
#> # ℹ 4 variables: taxon <chr>, method <chr>, matched_primers <chr>,
#> #   n_reads <dbl>
```
