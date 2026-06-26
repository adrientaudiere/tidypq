# Detect chimeric taxa with dada2 or vsearch

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Flags chimeric taxa using either
[`dada2::removeBimeraDenovo()`](https://rdrr.io/pkg/dada2/man/removeBimeraDenovo.html)
de novo detection (`method = "dada2"`) or UCHIME reference-based
detection via [vsearch](https://github.com/torognes/vsearch)
(`method = "vsearch_ref"`).

This is a pure detector: it flags taxa and returns a
[contam_tbl](https://adrientaudiere.github.io/tidypq/reference/contam_tbl.md)
without modifying `physeq`. Use
[`filter_contam_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_contam_pq.md)
to remove the flagged taxa.

## Usage

``` r
identify_contam_chimera_pq(
  physeq,
  method = c("dada2", "vsearch_ref"),
  dada2_method = "consensus",
  database = NULL,
  vsearchpath = "vsearch",
  vsearch_args = "--fasta_width 0",
  keep_temporary_files = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object with a populated `@refseq` slot
  containing DNA sequences.

- method:

  (character, default `"dada2"`) Detection engine: `"dada2"` for de novo
  abundance-based detection, or `"vsearch_ref"` for reference-based
  UCHIME detection (requires `database`).

- dada2_method:

  (character, default `"consensus"`) Passed as `method` to
  [`dada2::removeBimeraDenovo()`](https://rdrr.io/pkg/dada2/man/removeBimeraDenovo.html)
  when `method = "dada2"`. One of `"consensus"`, `"pooled"`,
  `"per-sample"`.

- database:

  (character, default `NULL`) Path to the reference FASTA database used
  by vsearch. Required when `method = "vsearch_ref"`.

- vsearchpath:

  (character, default `"vsearch"`) Path to the vsearch executable.

- vsearch_args:

  (character, default `"--fasta_width 0"`) Extra arguments passed
  verbatim to vsearch.

- keep_temporary_files:

  (logical, default `FALSE`) Keep the temporary FASTA files written to
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html) (vsearch only).

- verbose:

  (logical, default `TRUE`) Print a one-line detection summary.

- ...:

  Additional arguments passed to
  [`dada2::removeBimeraDenovo()`](https://rdrr.io/pkg/dada2/man/removeBimeraDenovo.html).

## Value

A
[contam_tbl](https://adrientaudiere.github.io/tidypq/reference/contam_tbl.md)
with one row per chimeric taxon, `method = "chimera"`, and the evidence
columns `chimera_method` (`"dada2"` or `"vsearch_ref"`) and `n_reads`.
Empty when no chimera is detected.

## See also

[`filter_contam_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_contam_pq.md)
to remove flagged taxa,
[`create_chimera_pq()`](https://adrientaudiere.github.io/tidypq/reference/create_chimera_pq.md)
to generate synthetic chimeras for benchmarking.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
data(data_fungi)
# \donttest{
flagged <- identify_contam_chimera_pq(data_fungi, method = "dada2")
#> ℹ chimera (dada2): flagged 450 of 1420 taxa.
flagged
#> contam_tbl: 450 flagged taxons from methods "chimera"
#> # A tibble: 450 × 4
#>    taxon  method  chimera_method n_reads
#>    <chr>  <chr>   <chr>            <dbl>
#>  1 ASV94  chimera dada2             5842
#>  2 ASV120 chimera dada2             4558
#>  3 ASV128 chimera dada2             4298
#>  4 ASV158 chimera dada2             3193
#>  5 ASV159 chimera dada2             3164
#>  6 ASV170 chimera dada2             2860
#>  7 ASV178 chimera dada2             2743
#>  8 ASV182 chimera dada2             2712
#>  9 ASV186 chimera dada2             2611
#> 10 ASV201 chimera dada2             2294
#> # ℹ 440 more rows
# }
```
