# Filter phyloseq by tree topology

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

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
# \donttest{
library(MiscMetabar)
library(ape)
#> 
#> Attaching package: ‘ape’
#> The following object is masked from ‘package:dplyr’:
#> 
#>     where
data(data_fungi)
pg <- data_fungi
pg@phy_tree <- ape::rtree(phyloseq::ntaxa(data_fungi),
                           tip.label = phyloseq::taxa_names(data_fungi))
filter_tree_pq(pg, taxa = phyloseq::taxa_names(data_fungi)[1:10])
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 10 taxa and 155 samples ]
#> sample_data() Sample Data:       [ 155 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 10 taxa by 12 taxonomic ranks ]
#> phy_tree()    Phylogenetic Tree: [ 10 tips and 9 internal nodes ]
#> refseq()      DNAStringSet:      [ 10 reference sequences ]
# }
```
