# Create a phyloseq object with synthetic chimeric sequences

This function creates synthetic chimeric sequences by combining parts of
existing sequences from a phyloseq object. Useful for benchmarking
chimera detection methods like
[`MiscMetabar::chimera_removal_vs()`](https://adrientaudiere.github.io/MiscMetabar/reference/chimera_removal_vs.html)
or
[`chimera_removal_dada2()`](https://adrientaudiere.github.io/tidypq/reference/chimera_removal_dada2.md).

## Usage

``` r
create_chimera_pq(
  physeq,
  n_chimeras = 5,
  prop_mean = 0.5,
  prop_sd = 0.15,
  prop_min = 0.1,
  seed = 123,
  median_abundance_multiplier = 0.1,
  ensure_distinct_parents = TRUE,
  min_parent_distance = 0.1
)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object with a refseq slot containing
  DNA sequences.

- n_chimeras:

  (integer, default: 5) Number of chimeric sequences to create.

- prop_mean:

  (numeric, default: 0.5) Mean of the normal distribution used to sample
  the proportion of the first parent sequence. A value of 0.5 means
  chimeras will be centered around 50/50 splits.

- prop_sd:

  (numeric, default: 0.15) Standard deviation of the normal distribution
  used to sample proportions. Higher values create more variable chimera
  breakpoints.

- prop_min:

  (numeric, default: 0.1) Minimum proportion threshold. Proportions
  below this value (or above 1 - prop_min) are resampled to ensure each
  parent contributes meaningfully to the chimera.

- seed:

  (integer, default: 123) Random seed for reproducibility.

- median_abundance_multiplier:

  (numeric, default: 0.1) Multiplier to set the abundance of chimeric
  sequences relative to the median abundance of existing sequences. A
  value of 0.1 means chimeras will have approximately 10% of the median
  abundance.

- ensure_distinct_parents:

  (logical, default: TRUE) If TRUE, ensures that parent2 is sufficiently
  different from parent1 based on `min_parent_distance`. Chimeras
  created from very similar parent sequences may be undetectable.

- min_parent_distance:

  (numeric, default: 0.1) Minimum sequence distance (proportion of
  differing positions) between parent1 and parent2. Only used when
  `ensure_distinct_parents = TRUE`. Values typically range from 0.05 (5%
  divergence) to 0.3 (30% divergence).

## Value

A list containing:

- physeq:

  The new phyloseq object with added chimeric sequences

- chimera_names:

  Character vector of chimera taxa names

- parent_info:

  Data frame with details about each chimera: chimera name, parent1,
  parent2, parent_distance, prop_parent1, breakpoint, seq_length

- params:

  List of parameters used (prop_mean, prop_sd, prop_min,
  ensure_distinct_parents, min_parent_distance)

## See also

[`MiscMetabar::chimera_removal_vs()`](https://adrientaudiere.github.io/MiscMetabar/reference/chimera_removal_vs.html),
[`chimera_removal_dada2()`](https://adrientaudiere.github.io/tidypq/reference/chimera_removal_dada2.md)

## Author

Adrien Taudiere

## Examples

``` r
if (FALSE) { # \dontrun{
library(MiscMetabar)
data(data_fungi)

# Default: centered around 50% with some variation
result <- create_chimera_pq(data_fungi, n_chimeras = 40)
data_fungi_test <- result$physeq
known_chimeras <- result$chimera_names

# View the parent information and proportions
print(result$parent_info)

# More variable proportions (wider distribution)
result2 <- create_chimera_pq(data_fungi, n_chimeras = 40,
                             prop_mean = 0.5, prop_sd = 0.25)

# Biased toward more of parent1 (e.g., 70/30 splits on average)
result3 <- create_chimera_pq(data_fungi, n_chimeras = 40,
                             prop_mean = 0.7, prop_sd = 0.1)

# Benchmark chimera detection methods
if (MiscMetabar::is_vsearch_installed()) {
  nochim_vs <- MiscMetabar::chimera_removal_vs(data_fungi_test)
  detected_vs <- known_chimeras[!known_chimeras %in% phyloseq::taxa_names(nochim_vs)]
  cat("vsearch detected:", length(detected_vs), "/",
      length(known_chimeras), "chimeras\n")
}

# Visualize the distribution of proportions
hist(result$parent_info$prop_parent1,
     main = "Distribution of parent1 proportions",
     xlab = "Proportion from parent1", xlim = c(0, 1))

# Ensure parents are at least 15% different (more detectable chimeras)
result4 <- create_chimera_pq(data_fungi, n_chimeras = 40,
                             min_parent_distance = 0.15)

# Disable parent distance filtering (allows similar parents)
result5 <- create_chimera_pq(data_fungi, n_chimeras = 40,
                             ensure_distinct_parents = FALSE)
} # }
```
