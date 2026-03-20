# Rename columns in tax_table

Rename columns (taxonomic ranks) in tax_table using tidyselect
semantics.

## Usage

``` r
rename_taxa_pq(physeq, ...)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- ...:

  Name-value pairs where the name is the new name and the value is the
  old name. Use `new_name = old_name` syntax.

## Value

A phyloseq object with renamed tax_table columns.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
# Rename a single rank
rename_taxa_pq(data_fungi, tax_domain = Domain)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]

# Rename multiple ranks
rename_taxa_pq(data_fungi, tax_phylum = Phylum, tax_class = Class)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]
```
