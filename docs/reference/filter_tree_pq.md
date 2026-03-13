# Filter phyloseq by tree topology

Filter a phyloseq object to include only taxa that are present in the
phylogenetic tree, or prune the tree to match the taxa in the phyloseq.
Can also filter based on tree properties like tip labels matching a
pattern.

## Usage

``` r
filter_tree_pq(
  physeq,
  taxa = NULL,
  pattern = NULL,
  invert = FALSE,
  clean_phyloseq_object = TRUE
)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object with a phy_tree slot.

- taxa:

  (character, optional) Character vector of taxa names to keep. If NULL,
  keeps all taxa that are present in both the OTU table and tree.

- pattern:

  (character, optional) A regular expression pattern to match against
  tip labels. Only tips matching the pattern are kept.

- invert:

  (logical, default FALSE) If TRUE and pattern is provided, keep tips
  that do NOT match the pattern.

- clean_phyloseq_object:

  if TRUE (default), the resulting phyloseq object is cleaned using
  [`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.html)
  to remove empty taxa/samples.

## Value

A phyloseq object with filtered tree and matching taxa.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
# Synchronize tree with OTU table (remove taxa not in tree)
# filter_tree_pq(physeq_with_tree)

# Keep only specific taxa in tree
# filter_tree_pq(physeq_with_tree, taxa = c("ASV1", "ASV2", "ASV3"))

# Filter by tip label pattern
# filter_tree_pq(physeq_with_tree, pattern = "^ASV")
```
