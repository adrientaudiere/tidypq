# tidypq

Tidyverse-style verbs for manipulating phyloseq objects at four scales: **samples**, **taxa**, **occurrences**, and **tree**.

## Function Naming Convention

All functions follow the pattern: `{verb}_{scale}_pq()`

| Verb | samples | taxa | occurrences | tree |
|------|---------|------|-------------|------|
| **filter** | `filter_samples_pq` | `filter_taxa_pq` | `filter_occurrences_pq` | `filter_tree_pq` |
| **select** | `select_samdata_pq` | `select_taxa_pq` | - | - |
| **mutate** | `mutate_samdata_pq` | `mutate_taxa_pq` | `mutate_occurrences_pq` | - |
| **slice** | `slice_samples_pq` | `slice_taxa_pq` | - | - |
| **arrange** | `arrange_samples_pq` | `arrange_taxa_pq` | - | - |
| **rename** | `rename_samples_pq` | `rename_taxa_pq` | - | - |

**Helper:** `taxa_prevalence(physeq, threshold)`

## Installation

```r
# Install from CRAN
install.packages("tidypq")

# Or install the development version from GitHub
devtools::install_github("adrientaudiere/tidypq")
```

## Usage

```r
library(tidypq)
library(MiscMetabar)

# Use the `.` pronoun to refer to the phyloseq object
data_fungi |>
filter_samples_pq(Height == "Low", sample_sums(.) > 1000) |>
filter_taxa_pq(taxa_sums(.) > 100) |>
mutate_samdata_pq(log_depth = log(sample_sums(.))) |>
mutate_occurrences_pq(. / sample_total)  # Relative abundance
```
