# Getting Started with tidypq

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

## Overview

tidypq provides tidyverse-style verbs for manipulating phyloseq objects.
Use the `.` pronoun to refer to the phyloseq object within expressions.

## Filter samples

``` r
# By metadata
data_fungi |> 
  filter_samples_pq(Height == "Low")
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 963 taxa and 45 samples ]
#> sample_data() Sample Data:       [ 45 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 963 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 963 reference sequences ]

# By sequencing depth
data_fungi |> 
  filter_samples_pq(sample_sums(.) > 5000)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1383 taxa and 108 samples ]
#> sample_data() Sample Data:       [ 108 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1383 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1383 reference sequences ]

# Both criteria
data_fungi |> 
  filter_samples_pq(Height == "Low" & sample_sums(.) > 5000)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 816 taxa and 20 samples ]
#> sample_data() Sample Data:       [ 20 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 816 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 816 reference sequences ]
```

## Filter taxa

``` r
# By taxonomy
data_fungi |> 
  filter_taxa_pq(Phylum == "Basidiomycota")
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 345 taxa and 181 samples ]
#> sample_data() Sample Data:       [ 181 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 345 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 345 reference sequences ]

# By abundance
data_fungi |> 
  filter_taxa_pq(taxa_sums(.) > 100)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1078 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1078 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1078 reference sequences ]
```

## Mutate

``` r
# Add sample depth
fungi_modif <- data_fungi |> 
  mutate_samdata_pq(depth = sample_sums(.), 
    hill_2 = vegan::renyi(.@otu_table, scales = 2, hill = TRUE),
    )
```

## Select columns

``` r
# Select sample metadata
data_fungi |> select_samdata_pq(Height, Time)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 2 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]

# Select taxonomic ranks
data_fungi |> select_taxa_pq(Phylum:Genus)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 5 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]
```

## Chaining operations

``` r
data_fungi |>
  filter_samples_pq(Height == "Low") |>
  filter_taxa_pq(taxa_sums(.) > 50) |>
  mutate_samdata_pq(log_depth = log(sample_sums(.)))
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 419 taxa and 45 samples ]
#> sample_data() Sample Data:       [ 45 samples by 8 sample variables ]
#> tax_table()   Taxonomy Table:    [ 419 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 419 reference sequences ]
```
