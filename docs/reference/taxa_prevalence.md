# Calculate taxa prevalence

Calculate the prevalence (number of samples in which a taxon is present)
for each taxon in a phyloseq object.

## Usage

``` r
taxa_prevalence(physeq, threshold = 0)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- threshold:

  (numeric, default 0) Minimum abundance to consider a taxon as present
  in a sample. Values \> threshold are counted as present.

## Value

A named numeric vector (class integer) with prevalence for each taxon.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
# Get prevalence for each taxon
prev <- taxa_prevalence(data_fungi)
#> Taxa are now in rows.
head(prev)
#>  ASV2  ASV6  ASV7  ASV8 ASV10 ASV12 
#>   137    25    44    92    48    34 

# Get prevalence with minimum abundance of 10
prev10 <- taxa_prevalence(data_fungi, threshold = 10)
#> Taxa are now in rows.
head(prev10)
#>  ASV2  ASV6  ASV7  ASV8 ASV10 ASV12 
#>    35    11    36    35    33    14 

# Use in filter_taxa_pq
# Keep taxa present in at least 5 samples
filter_taxa_pq(data_fungi, taxa_prevalence(.) >= 5)
#> Taxa are now in rows.
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 817 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 817 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 817 reference sequences ]
```
