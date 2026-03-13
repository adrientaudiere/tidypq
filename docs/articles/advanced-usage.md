# Advanced Usage

``` r
library(tidypq)
#> Loading required package: phyloseq
#> Loading required package: MiscMetabar
#> Loading required package: ggplot2
#> Loading required package: dada2
#> Loading required package: Rcpp
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
#> Loading required package: purrr
library(MiscMetabar)
```

## Custom functions in filter

``` r
# Filter taxa present in at least 3 samples
data_fungi |> 
  filter_taxa_pq(taxa_prevalence(.) >= 3)
#> Taxa are now in rows.
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1064 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1064 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1064 reference sequences ]

# Filter samples with above-median depth
data_fungi |>
   filter_samples_pq(sample_sums(.) > median(sample_sums(.)))
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1353 taxa and 92 samples ]
#> sample_data() Sample Data:       [ 92 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1353 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1353 reference sequences ]
```

## Custom functions in mutate

``` r
# Prevalence ratio
data_fungi |> mutate_taxa_pq(
  prev_ratio = taxa_prevalence(.) / nsamples(.)
)
#> Taxa are now in rows.
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 13 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]

# Relative abundance rank
data_fungi |> mutate_taxa_pq(
  abundance_rank = rank(-taxa_sums(.))
)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 13 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]
```

## Occurrence-level transformations

``` r
# Relative abundance
data_fungi |> mutate_occurrences_pq(. / sample_total)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]

# Filter low-abundance occurrences
data_fungi |> filter_occurrences_pq(. / sample_total > 0.001)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1387 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1387 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1387 reference sequences ]
```

## Complex pipelines

``` r
data_fungi |>
  filter_samples_pq(sample_sums(.) > quantile(sample_sums(.), 0.25)) |>
  filter_taxa_pq(taxa_prevalence(.) >= 2, taxa_sums(.) > 10) |>
  mutate_taxa_pq(
    prevalence = taxa_prevalence(.),
    mean_abundance = taxa_sums(.) / nsamples(.)
  ) |>
  mutate_occurrences_pq(. / sample_total)
#> Taxa are now in rows.
#> Taxa are now in rows.
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1185 taxa and 138 samples ]
#> sample_data() Sample Data:       [ 138 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1185 taxa by 14 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1185 reference sequences ]
```
