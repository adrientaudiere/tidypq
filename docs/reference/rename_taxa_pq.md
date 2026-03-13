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
rename_taxa_pq(data_fungi, tax_kingdom = Kingdom)
#> Error in dplyr::rename(tax_df, ...): Can't rename columns that don't exist.
#> ✖ Column `Kingdom` doesn't exist.

# Rename multiple ranks
rename_taxa_pq(data_fungi, tax_kingdom = Kingdom, tax_phylum = Phylum)
#> Error in dplyr::rename(tax_df, ...): Can't rename columns that don't exist.
#> ✖ Column `Kingdom` doesn't exist.
```
