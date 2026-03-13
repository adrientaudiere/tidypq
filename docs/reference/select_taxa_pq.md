# Select columns from tax_table in a phyloseq object

Select tax_table columns (taxonomic ranks) using tidyselect semantics.

## Usage

``` r
select_taxa_pq(physeq, ...)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- ...:

  One or more unquoted expressions separated by commas. Variable names
  can be used as if they were positions in the data frame, so
  expressions like `Kingdom:Genus` can be used to select a range of
  taxonomic ranks.

## Value

A phyloseq object with selected tax_table columns.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
# Select specific ranks
select_taxa_pq(data_fungi, Kingdom, Phylum, Class)
#> Error in select_taxa_pq(data_fungi, Kingdom, Phylum, Class): Can't select columns that don't exist.
#> ✖ Column `Kingdom` doesn't exist.

# Select a range of ranks
select_taxa_pq(data_fungi, Kingdom:Genus)
#> Error in select_taxa_pq(data_fungi, Kingdom:Genus): Can't select columns that don't exist.
#> ✖ Column `Kingdom` doesn't exist.

# Exclude ranks
select_taxa_pq(data_fungi, !Species)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 11 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]
```
