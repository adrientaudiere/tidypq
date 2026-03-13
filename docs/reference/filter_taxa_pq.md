# Filter taxa in a phyloseq object

Filter taxa using data masking on tax_table. Supports the `.` pronoun to
refer to the phyloseq object for use with functions like
[`taxa_sums()`](https://rdrr.io/pkg/phyloseq/man/taxa_sums.html).

## Usage

``` r
filter_taxa_pq(physeq, ..., clean_phyloseq_object = TRUE)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- ...:

  Expressions that return a logical value, evaluated in the context of
  tax_table. Multiple conditions are combined with `&`. Use `.` to refer
  to the phyloseq object (e.g., `taxa_sums(.) > 100`).

- clean_phyloseq_object:

  if TRUE (default), the resulting phyloseq object is cleaned using
  [`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.html)
  to remove empty taxa/samples.

## Value

A phyloseq object with filtered taxa.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
# Filter by taxonomy
filter_taxa_pq(data_fungi, Phylum == "Basidiomycota")
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 345 taxa and 181 samples ]
#> sample_data() Sample Data:       [ 181 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 345 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 345 reference sequences ]

# Filter by total abundance
filter_taxa_pq(data_fungi, taxa_sums(.) > 100)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1078 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1078 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1078 reference sequences ]

# Combine multiple conditions
filter_taxa_pq(data_fungi, Phylum == "Basidiomycota", taxa_sums(.) > 100)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 282 taxa and 180 samples ]
#> sample_data() Sample Data:       [ 180 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 282 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 282 reference sequences ]

# Keep taxa above median abundance
filter_taxa_pq(data_fungi, taxa_sums(.) > median(taxa_sums(.)))
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 709 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 709 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 709 reference sequences ]
```
