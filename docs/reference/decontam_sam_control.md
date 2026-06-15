# Decontaminate based on negative/blank control samples

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Remove potential contaminants by setting OTU values to 0 when they are
at or below the level observed in negative/blank control samples for
that particular OTU. If multiple control are available, for each taxon a
threshold is computed from the control samples using a summary function
(default: `max`). Occurrences in non-control samples that are at or
below this threshold are set to 0.

## Usage

``` r
decontam_sam_control(
  physeq,
  control_condition,
  fun = max,
  global_threshold = FALSE,
  remove_controls = FALSE,
  clean_phyloseq_object = TRUE,
  verbose = TRUE
)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- control_condition:

  An expression evaluated on sample_data that returns TRUE for control
  samples. Use `.` to refer to the phyloseq object (e.g.,
  `sample_type == "negative"`, `is_control == TRUE`).

- fun:

  (function, default `max`) A function to summarize the control sample
  values for each taxon (or globally if `global_threshold = TRUE`).
  Common choices: `max` (most conservative, default), `mean`, `median`,
  or a custom function. The function should take a numeric vector and
  return a single value.

- global_threshold:

  (logical, default FALSE) If TRUE, compute a single global threshold
  from all control occurrences instead of per-taxon thresholds. This
  applies `fun` to all values in the control samples.

- remove_controls:

  (logical, default FALSE) Whether to remove the control samples from
  the output phyloseq object after decontamination.

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
# Add a mock control column for demonstration using the 3 samples with lowest
# total abundance as controls
pq <- mutate_samdata_pq(data_fungi, is_control = sample_sums(.) < sort(sample_sums(.))[3])

# Decontaminate using max of controls as threshold (per-taxon)
decontam_sam_control(pq, is_control)
#> Decontamination complete.
#> Threshold type: per-taxon
#> Control samples: 2
#> Non-control samples: 183
#> Threshold function: max
#> Controls removed: FALSE
#> Sequences: 1839124 -> 1838970 (-154)
#> Occurrences: 12499 -> 12362 (-137)
#> Taxa: 1420 -> 1420 (-0)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 8 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]

# Use a global threshold (single value for all taxa)
decontam_sam_control(pq, is_control, global_threshold = TRUE)
#> Decontamination complete.
#> Threshold type: global
#> Control samples: 2
#> Non-control samples: 183
#> Threshold function: max
#> Global threshold value: 4
#> Controls removed: FALSE
#> Sequences: 1839124 -> 1828239 (-10885)
#> Occurrences: 12499 -> 6700 (-5799)
#> Taxa: 1420 -> 1410 (-10)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1410 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 8 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1410 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1410 reference sequences ]

# Use mean instead of max (less conservative)
decontam_sam_control(pq, is_control, fun = mean)
#> Decontamination complete.
#> Threshold type: per-taxon
#> Control samples: 2
#> Non-control samples: 183
#> Threshold function: mean
#> Controls removed: FALSE
#> Sequences: 1839124 -> 1839052 (-72)
#> Occurrences: 12499 -> 12430 (-69)
#> Taxa: 1420 -> 1420 (-0)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 8 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]

# Keep control samples in output
decontam_sam_control(pq, is_control, remove_controls = FALSE)
#> Decontamination complete.
#> Threshold type: per-taxon
#> Control samples: 2
#> Non-control samples: 183
#> Threshold function: max
#> Controls removed: FALSE
#> Sequences: 1839124 -> 1838970 (-154)
#> Occurrences: 12499 -> 12362 (-137)
#> Taxa: 1420 -> 1420 (-0)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 8 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]

# Use a custom function (e.g., 2x the max)
decontam_sam_control(pq, is_control, fun = \(x) 2 * max(x))
#> Decontamination complete.
#> Threshold type: per-taxon
#> Control samples: 2
#> Non-control samples: 183
#> Threshold function: function(x) 2 * max(x)
#> Controls removed: FALSE
#> Sequences: 1839124 -> 1838804 (-320)
#> Occurrences: 12499 -> 12294 (-205)
#> Taxa: 1420 -> 1420 (-0)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 8 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]
```
