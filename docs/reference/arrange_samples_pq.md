# Arrange samples by column values

Reorder samples based on sample_data columns. Supports the `.` pronoun
to refer to the phyloseq object for sorting by computed values.

## Usage

``` r
arrange_samples_pq(physeq, ..., clean_phyloseq_object = TRUE)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- ...:

  Variables to sort by. Use `desc()` for descending order. Use `.` to
  refer to the phyloseq object.

- clean_phyloseq_object:

  if TRUE (default), the resulting phyloseq object is cleaned using
  [`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.html)
  to remove empty taxa/samples.

## Value

A phyloseq object with reordered samples.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
# Arrange by a single column
arrange_samples_pq(data_fungi, Height)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 67 taxa and 1 samples ]
#> sample_data() Sample Data:       [ 1 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 67 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 67 reference sequences ]

# Arrange by sequencing depth (descending)
arrange_samples_pq(data_fungi, dplyr::desc(sample_sums(.)))
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 67 taxa and 1 samples ]
#> sample_data() Sample Data:       [ 1 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 67 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 67 reference sequences ]
```
