# Rename columns in sample_data

Rename columns in sample_data using tidyselect semantics.

## Usage

``` r
rename_samples_pq(physeq, ...)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- ...:

  Name-value pairs where the name is the new name and the value is the
  old name. Use `new_name = old_name` syntax.

## Value

A phyloseq object with renamed sample_data columns.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
# Rename a single column
rename_samples_pq(data_fungi, sample_height = Height)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]

# Rename multiple columns
rename_samples_pq(data_fungi, sample_height = Height, sample_time = Time)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]
```
