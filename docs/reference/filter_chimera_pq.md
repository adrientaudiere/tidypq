# Remove chimeric taxa in one step

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

A convenience wrapper that runs
[`identify_contam_chimera_pq()`](https://adrientaudiere.github.io/tidypq/reference/identify_contam_chimera_pq.md)
and then
[`filter_contam_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_contam_pq.md),
returning the chimera-free phyloseq object directly. Chimera removal is
a near-ubiquitous cleaning step, so it gets a dedicated verb; the other
contaminant methods are composed explicitly via their
`identify_contam_*_pq()` detector and
[`filter_contam_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_contam_pq.md).

## Usage

``` r
filter_chimera_pq(
  physeq,
  method = c("dada2", "vsearch_ref"),
  clean_phyloseq_object = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object with a populated `@refseq` slot
  containing DNA sequences.

- method:

  (character, default `"dada2"`) Detection engine, passed to
  [`identify_contam_chimera_pq()`](https://adrientaudiere.github.io/tidypq/reference/identify_contam_chimera_pq.md):
  `"dada2"` (de novo) or `"vsearch_ref"` (reference-based UCHIME,
  requires `database`).

- clean_phyloseq_object:

  (logical, default `TRUE`) Whether to clean the result with
  [`MiscMetabar::clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.html).

- verbose:

  (logical, default `TRUE`) Whether to report detection and removal.

- ...:

  Additional arguments passed to
  [`identify_contam_chimera_pq()`](https://adrientaudiere.github.io/tidypq/reference/identify_contam_chimera_pq.md)
  (e.g. `dada2_method`, `database`, `vsearchpath`, `vsearch_args`).

## Value

A phyloseq object with chimeric taxa removed (unchanged when none are
detected).

## See also

[`identify_contam_chimera_pq()`](https://adrientaudiere.github.io/tidypq/reference/identify_contam_chimera_pq.md)
to inspect chimeras before removal,
[`filter_contam_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_contam_pq.md)
for the general removal verb.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
data(data_fungi)
# \donttest{
data_nochim <- filter_chimera_pq(data_fungi, method = "dada2")
#> ℹ chimera (dada2): flagged 450 of 1420 taxa.
#> ! Removing 450 contaminant taxons flagged by methods "chimera".
# }
```
