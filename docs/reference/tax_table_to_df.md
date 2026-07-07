# Convert a `tax_table` into a data.frame with numeric and list columns

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

The `tax_table` slot of a phyloseq object is stored as a **character**
matrix, so numeric information (e.g. a bootstrap support, a trait value)
is silently coerced to text and list-valued cells cannot be represented
at all. This helper lifts those limitations by returning a tibble in
which columns that look numeric (or logical) are converted to their
natural type, and selected columns can be split into genuine list
columns.

It does **not** modify the phyloseq object; it produces a stand-alone
tibble for downstream (tidyverse) manipulation.

## Usage

``` r
tax_table_to_df(
  physeq,
  convert = TRUE,
  split = NULL,
  sep = "/",
  taxa_names_col = "taxon"
)
```

## Arguments

- physeq:

  (phyloseq or taxonomyTable, required) A phyloseq object with a
  `tax_table` slot, or a `tax_table` object directly.

- convert:

  (logical, default TRUE) If TRUE, each column is passed through
  [`utils::type.convert()`](https://rdrr.io/r/utils/type.convert.html)
  (with `as.is = TRUE`) so that columns holding only numbers become
  numeric/integer and `"TRUE"`/`"FALSE"` columns become logical. Columns
  with any non-numeric value are kept as character.

- split:

  (character, optional) Names of columns to turn into list columns by
  splitting each cell on `sep`. Useful when a rank cell packs several
  values (e.g. multiple GBIF matches).

- sep:

  (character, default `"/"`) Separator used to split the columns named
  in `split`. Passed to
  [`base::strsplit()`](https://rdrr.io/r/base/strsplit.html) with
  `fixed = TRUE`.

- taxa_names_col:

  (character or NULL, default `"taxon"`) Name of the leading column
  holding the taxa names. If NULL, taxa names are kept as row names
  instead of a column.

## Value

A [tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html)
with one row per taxon and one column per taxonomic rank, plus (unless
`taxa_names_col` is NULL) a leading column of taxa names.

## See also

[`phyloseq::tax_table()`](https://rdrr.io/pkg/phyloseq/man/tax_table-methods.html)

## Author

Adrien Taudière

## Examples

``` r
data_fungi_mini <-
 tidypq::mutate_taxa_pq(data_fungi_mini,
 Mol_Abundance = taxa_sums(.))

# Numeric-looking ranks are converted to numeric
df <- tax_table_to_df(data_fungi_mini)
df |>
  group_by(Class) |>
  summarize(mean=mean(Mol_Abundance))
#> # A tibble: 4 × 2
#>   Class                mean
#>   <chr>               <dbl>
#> 1 Agaricomycetes     13105.
#> 2 Atractiellomycetes 11053 
#> 3 Tremellomycetes     7750.
#> 4 NA                  5675 


if (FALSE) { # \dontrun{
# Split a packed column into a list column
  tax_table_to_df(data_fungi_mini, split = "Genus", sep = "_")
} # }
```
