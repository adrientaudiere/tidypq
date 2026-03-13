# Subset samples by position

Select samples by their integer positions.

## Usage

``` r
slice_samples_pq(physeq, ..., clean_phyloseq_object = TRUE)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- ...:

  Integer row indices. Positive values select samples, negative values
  drop samples.

- clean_phyloseq_object:

  if TRUE (default), the resulting phyloseq object is cleaned using
  [`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.html)
  to remove empty taxa/samples.

## Value

A phyloseq object with selected samples.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
# Select first 5 samples
slice_samples_pq(data_fungi, 1:5)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 268 taxa and 5 samples ]
#> sample_data() Sample Data:       [ 5 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 268 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 268 reference sequences ]

# Remove first 2 samples
slice_samples_pq(data_fungi, -(1:2))
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1415 taxa and 183 samples ]
#> sample_data() Sample Data:       [ 183 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1415 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1415 reference sequences ]
```
