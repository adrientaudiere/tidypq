# Add or modify columns in sample_data

Create new columns or modify existing ones in sample_data using data
masking. Supports the `.` pronoun to refer to the phyloseq object.

This function only modifies the sample_data slot (columns/metadata). It
cannot add or remove samples. The number of samples and sample names are
preserved.

## Usage

``` r
mutate_samdata_pq(physeq, ...)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- ...:

  Name-value pairs. The name gives the name of the column in the output.
  The value must be a vector of length 1 (recycled) or exactly the same
  length as the number of samples. Use `.` to refer to the

  phyloseq object.

## Value

A phyloseq object with modified sample_data (same samples, modified or
new columns).

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
# Add a new column based on sequencing depth
mutate_samdata_pq(data_fungi, log_depth = log(sample_sums(.)))
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 8 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]

# Modify an existing column
mutate_samdata_pq(data_fungi, Height = toupper(Height))
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]
```
