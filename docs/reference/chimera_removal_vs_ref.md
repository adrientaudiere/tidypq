# Remove chimeric sequences using vsearch and a reference database

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Applies UCHIME reference-based chimera detection via
[vsearch](https://github.com/torognes/vsearch) to identify and remove
chimeric sequences from a phyloseq object. Unlike de novo detection
([`MiscMetabar::chimera_removal_vs()`](https://adrientaudiere.github.io/MiscMetabar/reference/chimera_removal_vs.html)),
this method compares each query sequence against a trusted reference
database rather than inferring chimeras from abundance patterns within
the sample.

## Usage

``` r
chimera_removal_vs_ref(
  physeq,
  database,
  vsearchpath = "vsearch",
  vsearch_args = "--fasta_width 0",
  return_a_list = FALSE,
  keep_temporary_files = FALSE
)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object with a refseq slot containing
  DNA sequences.

- database:

  (character, required) Path to the reference FASTA database used by
  vsearch for chimera detection (e.g. a curated ITS or 16S database).

- vsearchpath:

  (character, default: `"vsearch"`) Path to the vsearch executable.
  Defaults to `"vsearch"` assuming it is on the system PATH.

- vsearch_args:

  (character, default: `"--fasta_width 0"`) Additional arguments passed
  verbatim to the vsearch command line.

- return_a_list:

  (logical, default: FALSE) If TRUE, returns a list with the filtered
  phyloseq, kept taxa names, and chimeric taxa names.

- keep_temporary_files:

  (logical, default: FALSE) If TRUE, temporary FASTA files written to
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html) are kept after the
  function returns.

## Value

If `return_a_list = FALSE` (default), returns a phyloseq object with
chimeric sequences removed. If `return_a_list = TRUE`, returns a list
with:

- physeq:

  The filtered phyloseq object

- kept_taxa:

  Character vector of retained taxa names

- chimeric_taxa:

  Character vector of removed chimeric taxa names

## Details

Sequences are written to a temporary FASTA file and passed to vsearch
using the `--uchime_ref` algorithm. Sequences classified as chimeras are
removed; sequences classified as non-chimeric or borderline are retained
(conservative default). Temporary files are cleaned up unless
`keep_temporary_files = TRUE`.

## See also

[`MiscMetabar::chimera_removal_vs()`](https://adrientaudiere.github.io/MiscMetabar/reference/chimera_removal_vs.html)
for de novo vsearch-based chimera removal,
[`chimera_removal_dada2()`](https://adrientaudiere.github.io/tidypq/reference/chimera_removal_dada2.md)
for dada2-based chimera removal,
[`create_chimera_pq()`](https://adrientaudiere.github.io/tidypq/reference/create_chimera_pq.md)
for creating test data with synthetic chimeras.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
data(data_fungi)

# Use the mini UNITE database bundled with MiscMetabar.
# NOTE: this is a very small subset of the full UNITE database,
# so a low number of chimeras detected here is expected and does not
# reflect what you would observe with the complete database.
mini_db <- system.file(
  "extdata", "mini_UNITE_fungi.fasta.gz",
  package = "MiscMetabar"
)

data_nochim <- chimera_removal_vs_ref(data_fungi, database = mini_db)
#> Discard 92 chimeric ASVs including:
#> ASV7;ASV26;ASV93;ASV94;ASV129;ASV170;ASV178;ASV209;ASV267;ASV328...
#> Try return_a_list=TRUE to see all discarded ASVs in the `chimeric_taxa` element.

# \donttest{
# Return detailed output including chimeric taxa names
result <- chimera_removal_vs_ref(
  data_fungi,
  database = mini_db,
  return_a_list = TRUE
)
#> Discard 92 chimeric ASVs including:
#> ASV7;ASV26;ASV93;ASV94;ASV129;ASV170;ASV178;ASV209;ASV267;ASV328...
#> Try return_a_list=TRUE to see all discarded ASVs in the `chimeric_taxa` element.
cat("Removed", length(result$chimeric_taxa), "chimeric ASVs\n")
#> Removed 92 chimeric ASVs

# Benchmark against synthetic chimeras created with create_chimera_pq()
pq_chim <- create_chimera_pq(data_fungi, n_chimeras = 10)
result_bench <- chimera_removal_vs_ref(
  pq_chim$physeq,
  database = mini_db,
  return_a_list = TRUE
)
#> Discard 100 chimeric ASVs including:
#> ASV7;ASV26;ASV93;ASV94;ASV129;ASV170;ASV178;ASV209;ASV267;ASV337...
#> Try return_a_list=TRUE to see all discarded ASVs in the `chimeric_taxa` element.
detected <- pq_chim$chimera_names[
  pq_chim$chimera_names %in% result_bench$chimeric_taxa
]
cat(
  "Detected", length(detected), "/", length(pq_chim$chimera_names),
  "synthetic chimeras\n"
)
#> Detected 8 / 10 synthetic chimeras
# }
```
