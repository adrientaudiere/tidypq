# Benchmarking Chimera Detection Methods

``` r
library(tidypq)
library(MiscMetabar)
library(phyloseq)
library(ggplot2)
library(patchwork)
```

## Overview

This vignette demonstrates how to benchmark chimera detection methods
using synthetic chimeric sequences. The
[`create_chimera_pq()`](https://adrientaudiere.github.io/tidypq/reference/create_chimera_pq.md)
function generates artificial chimeras with known parent sequences,
allowing you to evaluate detection sensitivity and specificity.

## Data Preparation

To reduce computation time, we use a subset of `data_fungi` filtered by
taxonomy. You can adjust the filtering to test on larger or smaller
datasets.

``` r
data(data_fungi)

# Filter to a single Class for faster computation
# Original data_fungi has 1420 taxa - this reduces it to ~200-400 taxa
data_fungi_subset <- subset_taxa(
  data_fungi,
  Class %in% c("Agaricomycetes", "Dacrymycetes","Eurotiomycetes", "Lecanoromycetes", "Tremellomycetes")
)

cat("Original taxa:", ntaxa(data_fungi), "\n")
cat("Subset taxa:", ntaxa(data_fungi_subset), "\n")
```

## Creating Test Data with Synthetic Chimeras

Create a phyloseq object with known chimeric sequences. The
`ensure_distinct_parents = TRUE` (default) ensures that parent sequences
are sufficiently different, making the chimeras detectable.

``` r
df_chim <- create_chimera_pq(
  data_fungi_subset,
  n_chimeras = 50,
  median_abundance_multiplier = 0.05,
  min_parent_distance = 0.05
)
data_fungi_test <- df_chim$physeq
known_chimeras <- df_chim$chimera_names

# View parent information including sequence distance
print(df_chim$parent_info)
cat(
  "\nMinimum parent distance:",
  min(df_chim$parent_info$parent_distance),
  "\n"
)
cat("\nMean parent distance:", mean(df_chim$parent_info$parent_distance), "\n")
```

## Comparing Detection Methods

### Helper function for benchmarking

``` r
# Helper function to benchmark a detection method
benchmark_detection <- function(
  physeq,
  method_fn,
  method_name,
  known_chimeras,
  ...
) {
  # Measure time and memory
  gc() # Clean up before measurement
  mem_before <- sum(gc()[, 2])

  start_time <- Sys.time()
  nochim <- method_fn(physeq, ...)
  end_time <- Sys.time()

  mem_after <- sum(gc()[, 2])

  # Calculate detection metrics
  detected <- known_chimeras[!known_chimeras %in% taxa_names(nochim)]
  missed <- known_chimeras[known_chimeras %in% taxa_names(nochim)]

  list(
    method = method_name,
    physeq_nochim = nochim,
    n_chimeras = length(known_chimeras),
    n_detected = length(detected),
    n_missed = length(missed),
    detection_rate = 100 * length(detected) / length(known_chimeras),
    total_removed = ntaxa(physeq) - ntaxa(nochim),
    time_seconds = as.numeric(difftime(end_time, start_time, units = "secs")),
    memory_mb = mem_after - mem_before,
    detected_names = detected,
    missed_names = missed
  )
}

# Helper function to calculate detection rate
calc_detection_rate <- function(physeq_nochim, known_chimeras) {
  detected <- sum(!known_chimeras %in% taxa_names(physeq_nochim))
  100 * detected / length(known_chimeras)
}
```

### Run benchmarks

``` r
# vsearch-based detection
result_vs <- benchmark_detection(
  data_fungi_test,
  chimera_removal_vs,
  "vsearch",
  known_chimeras
)

# dada2-based detection
result_dada2 <- benchmark_detection(
  data_fungi_test,
  chimera_removal_dada2,
  "dada2",
  known_chimeras
)
```

## Benchmark Summary Table

``` r
# Create summary table
benchmark_summary <- data.frame(
  Method = c(result_vs$method, result_dada2$method),
  `N Chimeras` = c(result_vs$n_chimeras, result_dada2$n_chimeras),
  Detected = c(result_vs$n_detected, result_dada2$n_detected),
  `Detection Rate (%)` = round(
    c(result_vs$detection_rate, result_dada2$detection_rate),
    1
  ),
  `Total Removed` = c(result_vs$total_removed, result_dada2$total_removed),
  `Time (s)` = round(c(result_vs$time_seconds, result_dada2$time_seconds), 2),
  `Memory (MB)` = round(c(result_vs$memory_mb, result_dada2$memory_mb), 1),
  check.names = FALSE
)

knitr::kable(benchmark_summary, caption = "Chimera Detection Benchmark Summary")
```

## Evaluating Detection Performance

``` r
cat("=== Chimera Detection Benchmark ===\n")
cat("Known chimeras:", length(known_chimeras), "\n\n")
cat(
  "vsearch detected:",
  result_vs$n_detected,
  "/",
  result_vs$n_chimeras,
  "(",
  round(result_vs$detection_rate, 1),
  "%)\n"
)
cat(
  "dada2 detected:",
  result_dada2$n_detected,
  "/",
  result_dada2$n_chimeras,
  "(",
  round(result_dada2$detection_rate, 1),
  "%)\n\n"
)

cat("vsearch total removed:", result_vs$total_removed, "\n")
cat("dada2 total removed:", result_dada2$total_removed, "\n")
```

## Analyzing Missed Chimeras

Examine which chimeras were missed and their parent sequences:

``` r
if (length(result_vs$missed_names) > 0) {
  cat("\nMissed by vsearch:\n")
  missed_info <- df_chim$parent_info[
    df_chim$parent_info$chimera %in% result_vs$missed_names,
  ]
  print(missed_info)
  cat(
    "\nMean parent distance of missed chimeras:",
    mean(missed_info$parent_distance),
    "\n"
  )
}

if (length(result_dada2$missed_names) > 0) {
  cat("\nMissed by dada2:\n")
  missed_info <- df_chim$parent_info[
    df_chim$parent_info$chimera %in% result_dada2$missed_names,
  ]
  print(missed_info)
  cat(
    "\nMean parent distance of missed chimeras:",
    mean(missed_info$parent_distance),
    "\n"
  )
}
```

## Comparing Detection Across Different Chimera Types

Create different types of synthetic chimeras to test detection
robustness:

``` r
# More variable proportions (harder to detect)
result_hard <- create_chimera_pq(
  data_fungi_subset,
  n_chimeras = 50,
  prop_mean = 0.5,
  prop_sd = 0.25 # wider distribution
)

# Biased proportions (70/30 splits)
result_biased <- create_chimera_pq(
  data_fungi_subset,
  n_chimeras = 50,
  prop_mean = 0.7,
  prop_sd = 0.1
)

# Lower abundance chimeras (harder to detect)
result_low_abund <- create_chimera_pq(
  data_fungi_subset,
  n_chimeras = 50,
  median_abundance_multiplier = 0.01
)
```

Run both detection methods on all chimera types:

``` r
# Run detection on all variants
datasets <- list(
  default = list(physeq = data_fungi_test, chimeras = known_chimeras),
  hard = list(
    physeq = result_hard$physeq,
    chimeras = result_hard$chimera_names
  ),
  biased = list(
    physeq = result_biased$physeq,
    chimeras = result_biased$chimera_names
  ),
  low_abund = list(
    physeq = result_low_abund$physeq,
    chimeras = result_low_abund$chimera_names
  )
)

results <- data.frame(
  dataset = character(),
  method = character(),
  detection_rate = numeric(),
  time_seconds = numeric(),
  stringsAsFactors = FALSE
)

for (name in names(datasets)) {
  ds <- datasets[[name]]

  # vsearch detection
  start_time <- Sys.time()
  nochim_vs <- chimera_removal_vs(ds$physeq)
  time_vs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  rate_vs <- calc_detection_rate(nochim_vs, ds$chimeras)

  # dada2 detection
  start_time <- Sys.time()
  nochim_dada2 <- chimera_removal_dada2(ds$physeq)
  time_dada2 <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  rate_dada2 <- calc_detection_rate(nochim_dada2, ds$chimeras)

  results <- rbind(
    results,
    data.frame(
      dataset = name,
      method = c("vsearch", "dada2"),
      detection_rate = c(rate_vs, rate_dada2),
      time_seconds = c(time_vs, time_dada2),
      stringsAsFactors = FALSE
    )
  )
}

print(results)
```

## Summary Table: All Scenarios

``` r
# Reshape for summary table
results_wide <- reshape(
  results,
  direction = "wide",
  idvar = "dataset",
  timevar = "method",
  v.names = c("detection_rate", "time_seconds")
)
names(results_wide) <- gsub("\\.", " ", names(results_wide))

knitr::kable(
  results_wide,
  caption = "Detection Rate and Time Across Scenarios",
  digits = 1
)
```

## Visualizing Results

``` r
# Combine parent info from all datasets for histograms
prop_data <- rbind(
  data.frame(
    dataset = "Default (sd=0.15)",
    prop = df_chim$parent_info$prop_parent1,
    min_dist = min(df_chim$parent_info$parent_distance)
  ),
  data.frame(
    dataset = "Hard (sd=0.25)",
    prop = result_hard$parent_info$prop_parent1,
    min_dist = min(result_hard$parent_info$parent_distance)
  ),
  data.frame(
    dataset = "Biased (mean=0.7)",
    prop = result_biased$parent_info$prop_parent1,
    min_dist = min(result_biased$parent_info$parent_distance)
  ),
  data.frame(
    dataset = "Low abundance",
    prop = result_low_abund$parent_info$prop_parent1,
    min_dist = min(result_low_abund$parent_info$parent_distance)
  )
)
prop_data$dataset <- factor(
  prop_data$dataset,
  levels = c(
    "Default (sd=0.15)",
    "Hard (sd=0.25)",
    "Biased (mean=0.7)",
    "Low abundance"
  )
)

# Distribution of parent1 proportions for each dataset
p_hist <- ggplot(prop_data, aes(x = prop)) +
  geom_histogram(
    binwidth = 0.05,
    fill = "steelblue",
    color = "white",
    alpha = 0.7
  ) +
  facet_wrap(~dataset, ncol = 2) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    x = "Proportion from parent1",
    y = "Count",
    title = "Distribution of Parent1 Proportions"
  ) +
  theme_bw()

print(p_hist)

# Barplot comparing all detection rates
results$dataset <- factor(
  results$dataset,
  levels = c(
    "Default (sd=0.15)",
    "Hard (sd=0.25)",
    "Biased (mean=0.7)",
    "Low abundance"
  )
)

p_bar <- ggplot(results, aes(x = dataset, y = detection_rate, fill = method)) +
  geom_violin(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(
    yintercept = seq(0, 100, 20),
    color = "gray",
    linetype = "dashed",
    alpha = 0.5
  ) +
  scale_fill_manual(values = c("vsearch" = "steelblue", "dada2" = "coral")) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(
    x = "Dataset",
    y = "Detection rate (%)",
    title = "Chimera Detection Across Different Scenarios",
    fill = "Method"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_bar)
```

## Effect of Chimera Abundance on Detection

Chimera abundance relative to parent sequences affects detection
sensitivity.

``` r
abundance_multipliers <- c(0.01, 0.05, 0.1, 0.5)

abundance_results <- data.frame(
  multiplier = numeric(),
  method = character(),
  detection_rate = numeric(),
  time_seconds = numeric(),
  stringsAsFactors = FALSE
)

for (mult in abundance_multipliers) {
  # Create chimeras with this abundance level
  result_abund <- create_chimera_pq(
    data_fungi_subset,
    n_chimeras = 50,
    median_abundance_multiplier = mult,
    seed = 123
  )

  # vsearch detection
  start_time <- Sys.time()
  nochim_vs <- chimera_removal_vs(result_abund$physeq)
  time_vs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  rate_vs <- calc_detection_rate(nochim_vs, result_abund$chimera_names)

  # dada2 detection
  start_time <- Sys.time()
  nochim_dada2 <- chimera_removal_dada2(result_abund$physeq)
  time_dada2 <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  rate_dada2 <- calc_detection_rate(nochim_dada2, result_abund$chimera_names)

  abundance_results <- rbind(
    abundance_results,
    data.frame(
      multiplier = c(mult, mult),
      method = c("vsearch", "dada2"),
      detection_rate = c(rate_vs, rate_dada2),
      time_seconds = c(time_vs, time_dada2),
      stringsAsFactors = FALSE
    )
  )
}

print(abundance_results)
```

### Abundance Summary Table

``` r
knitr::kable(
  abundance_results,
  caption = "Detection by Chimera Abundance Level",
  digits = c(2, 0, 1, 2)
)
```

### Visualizing Abundance Effect

``` r
# Detection rate vs abundance multiplier
p_abund_detection <- ggplot(
  abundance_results,
  aes(x = multiplier, y = detection_rate, color = method)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, 100)) +
  scale_color_manual(values = c("vsearch" = "steelblue", "dada2" = "coral")) +
  labs(
    x = "median_abundance_multiplier",
    y = "Detection rate (%)",
    title = "Detection vs Chimera Abundance",
    color = "Method"
  ) +
  theme_bw()

print(p_abund_detection)

# Barplot comparison
abundance_results$multiplier_label <- paste0(
  abundance_results$multiplier * 100,
  "%"
)
abundance_results$multiplier_label <- factor(
  abundance_results$multiplier_label,
  levels = paste0(sort(unique(abundance_results$multiplier)) * 100, "%")
)

p_abund_bar <- ggplot(
  abundance_results,
  aes(x = multiplier_label, y = detection_rate, fill = method)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(
    yintercept = seq(0, 100, 20),
    color = "gray",
    linetype = "dashed",
    alpha = 0.5
  ) +
  scale_fill_manual(values = c("vsearch" = "steelblue", "dada2" = "coral")) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(
    x = "Chimera abundance (% of median)",
    y = "Detection rate (%)",
    title = "Detection Rate by Chimera Abundance Level",
    fill = "Method"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_abund_bar)
```

## Effect of Parent Sequence Distance

Chimeras created from very similar parent sequences are harder to
detect. The `ensure_distinct_parents` parameter (default TRUE) helps
avoid this issue.

``` r
# Compare chimeras with different minimum parent distances
distance_results <- data.frame(
  min_param_distance = numeric(),
  method = character(),
  mean_parent_distance = numeric(),
  min_parent_distance = numeric(),
  max_parent_distance = numeric(),
  detection_rate = numeric(),
  nb_detected = numeric(),
  stringsAsFactors = FALSE
)

for (min_dist in c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.4)) {
  result_dist <- create_chimera_pq(
    data_fungi_subset,
    n_chimeras = 50,
    min_parent_distance = min_dist,
    seed = 456
  )

  mean_dist <- mean(result_dist$parent_info$parent_distance)
  min_dist_parent <- min(result_dist$parent_info$parent_distance)

  # vsearch detection
  nochim_vs <- chimera_removal_vs(result_dist$physeq)
  rate_vs <- calc_detection_rate(nochim_vs, result_dist$chimera_names)

  # dada2 detection
  nochim_dada2 <- chimera_removal_dada2(result_dist$physeq)
  rate_dada2 <- calc_detection_rate(nochim_dada2, result_dist$chimera_names)

  distance_results <- rbind(
    distance_results,
    data.frame(
      min_param_distance = min_dist,
      method = c("vsearch", "dada2"),
      mean_parent_distance = c(mean_dist, mean_dist),
      min_parent_distance = c(min_dist_parent, min_dist_parent),
      max_parent_distance = c(
        max(result_dist$parent_info$parent_distance),
        max(result_dist$parent_info$parent_distance)
      ),
      detection_rate = c(rate_vs, rate_dada2),
      nb_detected = c(
       ntaxa(result_dist$physeq) - ntaxa(nochim_vs),
        ntaxa(result_dist$physeq) - ntaxa(nochim_dada2)
      ),
      stringsAsFactors = FALSE
    )
  )
}

knitr::kable(
  distance_results,
  caption = "Detection Rate by Minimum Parent Distance",
  digits = c(2, 0, 3, 4, 0)
)
```

``` r
# Plot detection rate vs mean parent distance
p_dist <- ggplot(
  distance_results,
  aes(x = min_param_distance, y = detection_rate, color = method)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_continuous(limits = c(0, 0.5)) +
  scale_y_continuous(limits = c(0, 100)) +  
  scale_color_manual(values = c("vsearch" = "steelblue", "dada2" = "coral")) +
  labs(
    x = "Min parameter distance",
    y = "Detection rate (%)",
    title = "Detection vs Min Parameter  Distance",
    color = "Method"
  ) +
  theme_bw() +
  ggrepel::geom_label_repel(aes(label=nb_detected), size = 3, show.legend = FALSE)
print(p_dist)
```
