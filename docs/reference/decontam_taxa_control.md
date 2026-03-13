# Decontaminate based on control taxa

Remove potential contaminants by using known control taxa (e.g.,
spike-ins, synthetic sequences) to estimate background contamination
levels. For each sample, a threshold is computed from the control taxa
using a summary function (default: `max`). Occurrences of non-control
taxa that are at or below this threshold are set to 0.

## Usage

``` r
decontam_taxa_control(
  physeq,
  control_condition,
  fun = max,
  global_threshold = FALSE,
  remove_control_taxa = TRUE,
  clean_phyloseq_object = TRUE,
  verbose = TRUE
)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- control_condition:

  An expression evaluated on tax_table that returns TRUE for control
  taxa. Use `.` to refer to the phyloseq object. Examples:
  `Genus == "Tintelnotia"`, `Family == "Mitochondria"`,
  `taxa_names(.) %in% c("ASV1", "ASV2")`.

- fun:

  (function, default `max`) A function to summarize the control taxa
  values for each sample (or globally if `global_threshold = TRUE`).
  Common choices: `max` (most conservative, default), `mean`, `median`.

- global_threshold:

  (logical, default FALSE) If TRUE, compute a single global threshold
  from all control taxa occurrences instead of per-sample thresholds.

- remove_control_taxa:

  (logical, default TRUE) Whether to remove the control taxa from the
  output phyloseq object after decontamination.

- clean_phyloseq_object:

  (logical, default TRUE) Whether to clean the resulting phyloseq object
  using
  [`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.html)
  to remove empty taxa/samples.

- verbose:

  (logical, default TRUE) Whether to print additional information.

## Value

A phyloseq object with decontaminated OTU values.

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
# Using a condition on tax_table (e.g., select by Genus)
decontam_taxa_control(data_fungi, Genus == "Tintelnotia")
#> Decontamination complete.
#> Threshold type: per-sample
#> Number of control taxa: 1
#> Number of non-control taxa: 1419
#> Function used for threshold: max
#> Remove control taxa: TRUE
#> Number of sequences before: 1839124
#> Number of sequences after: 1835765
#> Number of sequences discarded: 3359
#> Number of occurrences discarded: 182
#> Number of taxa before: 1420
#> Number of taxa after: 1417
#> Number of discarded taxa: 3
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1417 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1417 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1417 reference sequences ]

# Using taxa names directly
control_taxa <- phyloseq::taxa_names(data_fungi)[1:2]
decontam_taxa_control(data_fungi, taxa_names(.) %in% control_taxa)
#> Decontamination complete.
#> Threshold type: per-sample
#> Number of control taxa: 2
#> Number of non-control taxa: 1418
#> Function used for threshold: max
#> Remove control taxa: TRUE
#> Number of sequences before: 1839124
#> Number of sequences after: 1578082
#> Number of sequences discarded: 261042
#> Number of occurrences discarded: 6193
#> Number of taxa before: 1420
#> Number of taxa after: 1351
#> Number of discarded taxa: 69
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1351 taxa and 169 samples ]
#> sample_data() Sample Data:       [ 169 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1351 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1351 reference sequences ]

# Use a global threshold
decontam_taxa_control(data_fungi, Genus == "Tintelnotia", global_threshold = TRUE)
#> Decontamination complete.
#> Threshold type: global
#> Number of control taxa: 1
#> Number of non-control taxa: 1419
#> Function used for threshold: max
#> Global threshold value: 210
#> Remove control taxa: TRUE
#> Number of sequences before: 1839124
#> Number of sequences after: 1615276
#> Number of sequences discarded: 223848
#> Number of occurrences discarded: 11470
#> Number of taxa before: 1420
#> Number of taxa after: 582
#> Number of discarded taxa: 838
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 582 taxa and 146 samples ]
#> sample_data() Sample Data:       [ 146 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 582 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 582 reference sequences ]

# Keep control taxa in output
decontam_taxa_control(data_fungi, Genus == "Tintelnotia", remove_control_taxa = FALSE)
#> Decontamination complete.
#> Threshold type: per-sample
#> Number of control taxa: 1
#> Number of non-control taxa: 1419
#> Function used for threshold: max
#> Remove control taxa: FALSE
#> Number of sequences before: 1839124
#> Number of sequences after: 1835977
#> Number of sequences discarded: 3147
#> Number of occurrences discarded: 180
#> Number of taxa before: 1420
#> Number of taxa after: 1418
#> Number of discarded taxa: 2
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1418 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1418 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1418 reference sequences ]
```
