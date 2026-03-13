# Filter occurrences in the OTU table

Set OTU table values to 0 based on a condition. This is useful for
removing singletons, low-abundance values, or other filtering operations
at the cell level.

The condition is evaluated vectorized across the entire OTU matrix with
access to special variables (all as matrices matching OTU dimensions):

- `.` = cell values (the OTU matrix)

- `sample_total` = sum of each sample (repeated per column)

- `taxon_total` = sum of each taxon (repeated per row)

- `sample_mean` = mean of each sample (repeated per column)

- `taxon_mean` = mean of each taxon (repeated per row)

Values that do not satisfy the condition are set to 0.

## Usage

``` r
filter_occurrences_pq(physeq, condition, clean_phyloseq_object = TRUE)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- condition:

  An expression evaluated on the OTU matrix. Values where the condition
  is FALSE (or NA) are set to 0. Use `.` to refer to cell values.

- clean_phyloseq_object:

  if TRUE (default), the resulting phyloseq object is cleaned using
  [`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.html)
  to remove empty taxa/samples.

## Value

A phyloseq object with filtered OTU values.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
# Remove singletons (abundance = 1)
filter_occurrences_pq(data_fungi, . > 1)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1416 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1416 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1416 reference sequences ]

# Keep only values above 0.01% of sample total
filter_occurrences_pq(data_fungi, . / sample_total > 0.0001)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1418 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1418 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1418 reference sequences ]

# Keep only values above taxon mean
filter_occurrences_pq(data_fungi, . > taxon_mean)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 184 samples ]
#> sample_data() Sample Data:       [ 184 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]
```
