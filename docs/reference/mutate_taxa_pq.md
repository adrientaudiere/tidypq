# Add or modify columns in tax_table

Create new columns or modify existing ones in tax_table using data
masking. Supports the `.` pronoun to refer to the phyloseq object.

This function only modifies the tax_table slot (columns/taxonomic
ranks). It cannot add or remove taxa. The number of taxa and taxa names
are preserved.

## Usage

``` r
mutate_taxa_pq(physeq, ...)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- ...:

  Name-value pairs. The name gives the name of the column in the output.
  The value must be a vector of length 1 (recycled) or exactly the same
  length as the number of taxa. Use `.` to refer to the phyloseq object.

## Value

A phyloseq object with modified tax_table (same taxa, modified or new
columns).

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
# Replace NA values in a column
mutate_taxa_pq(data_fungi, Genus = ifelse(is.na(Genus), "Unknown", Genus))
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]

# Add a new column based on abundance
mutate_taxa_pq(data_fungi, total_abundance = taxa_sums(.))
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 13 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]
```
