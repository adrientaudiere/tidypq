# Filter samples in a phyloseq object

Filter samples using data masking on sample_data. Supports the `.`
pronoun to refer to the phyloseq object for use with functions like
[`sample_sums()`](https://rdrr.io/pkg/phyloseq/man/sample_sums.html).

## Usage

``` r
filter_samples_pq(physeq, ..., clean_phyloseq_object = TRUE)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- ...:

  Expressions that return a logical value, evaluated in the context of
  sample_data. Multiple conditions are combined with `&`. Use `.` to
  refer to the phyloseq object (e.g., `sample_sums(.) > 1000`).

- clean_phyloseq_object:

  if TRUE (default), the resulting phyloseq object is cleaned using
  [`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.html)
  to remove empty taxa/samples.

## Value

A phyloseq object with filtered samples.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
# Filter by sample metadata
filter_samples_pq(data_fungi, Height == "Low")
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 963 taxa and 45 samples ]
#> sample_data() Sample Data:       [ 45 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 963 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 963 reference sequences ]

# Filter by sequencing depth
filter_samples_pq(data_fungi, sample_sums(.) > 1000)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1418 taxa and 143 samples ]
#> sample_data() Sample Data:       [ 143 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1418 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1418 reference sequences ]

# Combine multiple conditions
filter_samples_pq(data_fungi, Height == "Low", sample_sums(.) > 5000)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 816 taxa and 20 samples ]
#> sample_data() Sample Data:       [ 20 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 816 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 816 reference sequences ]

# Keep samples above median abundance
filter_samples_pq(data_fungi, sample_sums(.) > median(sample_sums(.)))
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1353 taxa and 92 samples ]
#> sample_data() Sample Data:       [ 92 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1353 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1353 reference sequences ]

# Keep samples above half of the average abundance
filter_samples_pq(data_fungi, sample_sums(.) > sum(sample_sums(.))/phyloseq::nsamples(.)/2)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1383 taxa and 108 samples ]
#> sample_data() Sample Data:       [ 108 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1383 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1383 reference sequences ]
```
