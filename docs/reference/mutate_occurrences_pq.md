# Transform OTU table values

Apply a transformation to all values in the OTU table. This is useful
for computing relative abundances, log transformations, or other
value-level operations.

The expression is evaluated vectorized across the entire OTU matrix with
access to special variables (all as matrices matching OTU dimensions):

- `.` = cell values (the OTU matrix)

- `sample_total` = sum of each sample (repeated per column)

- `taxon_total` = sum of each taxon (repeated per row)

- `sample_mean` = mean of each sample (repeated per column)

- `taxon_mean` = mean of each taxon (repeated per row)

- `sample_median` = median of each sample (repeated per column)

- `taxon_median` = median of each taxon (repeated per row)

## Usage

``` r
mutate_occurrences_pq(physeq, expr)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- expr:

  An expression evaluated on the OTU matrix. The result replaces the
  original values. Use `.` to refer to cell values.

## Value

A phyloseq object with transformed OTU values.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
# Convert to relative abundance (proportion of sample total)
mutate_occurrences_pq(data_fungi, . / sample_total)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]

# Log transformation (adding pseudocount)
mutate_occurrences_pq(data_fungi, log1p(.))
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]

# Center by taxon mean
mutate_occurrences_pq(data_fungi, . - taxon_mean)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]
```
