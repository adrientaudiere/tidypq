# Remove chimeric sequences using dada2

Applies
[`dada2::removeBimeraDenovo()`](https://rdrr.io/pkg/dada2/man/removeBimeraDenovo.html)
to identify and remove chimeric sequences from a phyloseq object based
on sequence abundance patterns.

## Usage

``` r
chimera_removal_dada2(physeq, method = "consensus", return_a_list = FALSE, ...)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object with a refseq slot containing
  DNA sequences.

- method:

  (character, default: "consensus") Method for chimera detection. Passed
  to
  [`dada2::removeBimeraDenovo()`](https://rdrr.io/pkg/dada2/man/removeBimeraDenovo.html).
  Options: "consensus", "pooled", "per-sample".

- return_a_list:

  (logical, default: FALSE) If TRUE, returns a list with the filtered
  phyloseq, kept taxa names, and chimeric taxa names.

- ...:

  Additional arguments passed to
  [`dada2::removeBimeraDenovo()`](https://rdrr.io/pkg/dada2/man/removeBimeraDenovo.html).

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

This function extracts sequences and their abundances from the phyloseq
object, applies dada2's de novo chimera detection algorithm, and returns
a pruned phyloseq object containing only non-chimeric sequences.

The dada2 method uses sequence abundance information to identify
chimeras, assuming that chimeric sequences are less abundant than their
parent sequences.

## See also

[`MiscMetabar::chimera_removal_vs()`](https://adrientaudiere.github.io/MiscMetabar/reference/chimera_removal_vs.html)
for vsearch-based chimera removal,
[`create_chimera_pq()`](https://adrientaudiere.github.io/tidypq/reference/create_chimera_pq.md)
for creating test data with synthetic chimeras.

## Examples

``` r
if (FALSE) { # \dontrun{
library(MiscMetabar)
data(data_fungi)

# Basic usage
data_nochim <- chimera_removal_dada2(data_fungi)

# Get detailed output
result <- chimera_removal_dada2(data_fungi, return_a_list = TRUE)
cat("Removed", length(result$chimeric_taxa), "chimeric ASVs\n")

# Use pooled method
data_nochim <- chimera_removal_dada2(data_fungi, method = "pooled")
} # }
```
