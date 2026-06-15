# Remove negative-control contaminants from a phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Drop every taxon that is detected above a read threshold in negative
controls, then optionally drop the negative-control samples themselves.

This is a simpler, more aggressive alternative to
[`decontam_sam_control()`](https://adrientaudiere.github.io/tidypq/reference/decontam_sam_control.md):
instead of zeroing OTU values at or below the per-taxon control
threshold, it removes the offending taxa altogether. Use it when any
presence in a negative control is considered disqualifying.

Algorithm:

1.  Subset `physeq` to its negative-control samples.

2.  Identify taxa whose summed reads across those samples is **strictly
    greater than** `min_reads`.

3.  Remove those taxa from the full dataset.

4.  If `drop_neg_samples = TRUE` (default), also drop the
    negative-control samples.

## Usage

``` r
neg_control_clean_pq(
  physeq,
  neg_control,
  min_reads = 2,
  drop_neg_samples = TRUE,
  clean_phyloseq_object = TRUE,
  verbose = TRUE
)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- neg_control:

  An expression evaluated on sample_data that returns TRUE for
  negative-control samples (e.g., `is_control == TRUE`,
  `sample_type == "blank"`). Use `.` to refer to the phyloseq object.

- min_reads:

  (integer, default: 2) Read-count threshold for a taxon to count as a
  contaminant. A taxon is removed when its total reads across NC samples
  is strictly greater than `min_reads`. The default `2` matches the
  common rule "tolerate up to two stray reads in NCs".

- drop_neg_samples:

  (logical, default: TRUE) Drop the negative-control samples from the
  returned object. Set to FALSE to keep them (e.g., for NC-aware
  downstream modelling).

- clean_phyloseq_object:

  (logical, default: TRUE) Whether to clean the resulting phyloseq
  object using
  [`MiscMetabar::clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.html)
  to remove empty taxa/samples.

- verbose:

  (logical, default: TRUE) Whether to print a short message reporting
  the number of contaminant taxa removed.

## Value

A phyloseq object with contaminant taxa (and, by default, NC samples)
removed.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
data(data_fungi)
# Mark the three lowest-depth samples as mock controls for demo
pq <- mutate_samdata_pq(
  data_fungi,
  is_control = sample_sums(.) < sort(sample_sums(.))[4]
)
neg_control_clean_pq(pq, is_control)
#> Removing 2 contaminant taxa (>2 reads in NCs).
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1418 taxa and 182 samples ]
#> sample_data() Sample Data:       [ 182 samples by 8 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1418 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1418 reference sequences ]

# Keep NC samples in output
neg_control_clean_pq(pq, is_control, drop_neg_samples = FALSE)
#> Removing 2 contaminant taxa (>2 reads in NCs).
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1418 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 8 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1418 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1418 reference sequences ]

# Stricter: drop any taxon ever detected in NCs
neg_control_clean_pq(pq, is_control, min_reads = 0)
#> Removing 14 contaminant taxa (>0 reads in NCs).
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1406 taxa and 182 samples ]
#> sample_data() Sample Data:       [ 182 samples by 8 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1406 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1406 reference sequences ]
```
