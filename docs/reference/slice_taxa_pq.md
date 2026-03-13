# Subset taxa by position

Select taxa by their integer positions.

## Usage

``` r
slice_taxa_pq(physeq, ..., clean_phyloseq_object = TRUE)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- ...:

  Integer row indices. Positive values select taxa, negative values drop
  taxa.

- clean_phyloseq_object:

  if TRUE (default), the resulting phyloseq object is cleaned using
  [`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.html)
  to remove empty taxa/samples.

## Value

A phyloseq object with selected taxa.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
# Select first 10 taxa
slice_taxa_pq(data_fungi, 1:10)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 10 taxa and 155 samples ]
#> sample_data() Sample Data:       [ 155 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 10 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 10 reference sequences ]

# Remove first 5 taxa
slice_taxa_pq(data_fungi, -(1:5))
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1415 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1415 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1415 reference sequences ]
```
