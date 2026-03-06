# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

**tidypq** is an R package providing tidyverse-style verbs for manipulating phyloseq objects at four scales: samples, taxa, occurrences, and tree. All functions follow the naming convention `{verb}_{scale}_pq()`.

## Common Commands

```bash
# Run code with loaded package
Rscript -e "devtools::load_all(); code"

# Run all tests
Rscript -e "devtools::test()"

# Run tests for files starting with {name}
Rscript -e "devtools::test(filter = '^{name}')"

# Run tests for R/{name}.R
Rscript -e "devtools::test_active_file('R/{name}.R')"

# Generate documentation
Rscript -e "devtools::document()"

# Full package check
Rscript -e "devtools::check()"
```

## Architecture

### Core Dependencies
- **phyloseq**: S4 class system for microbiome data (slots: @tax_table, @otu_table, @refseq, @sam_data)
- **MiscMetabar**: Parent package with metabarcoding utilities (verify_pq, clean_pq)
- **rlang**: Data masking and tidy evaluation

### Key Modules

| Module | Purpose |
|--------|---------|
| `tidypq-package.R` | Package documentation and imports |
| `utils.R` | Internal helpers for data masking |
| `samples.R` | filter, select, mutate, slice, arrange, rename for samples |
| `taxa.R` | filter, select, mutate, slice, arrange, rename for taxa |
| `occurrences.R` | filter, mutate for OTU table values |
| `tree.R` | filter for phylogenetic tree |
| `helpers.R` | taxa_prevalence and other helpers |

### Function Naming Convention

All functions follow: `{verb}_{scale}_pq()`

| Scale | filter | select | mutate | slice | arrange | rename |
|-------|--------|--------|--------|-------|---------|--------|
| samples | `filter_samples_pq` | `select_samdata_pq` | `mutate_samdata_pq` | `slice_samples_pq` | `arrange_samples_pq` | `rename_samples_pq` |
| taxa | `filter_taxa_pq` | `select_taxa_pq` | `mutate_taxa_pq` | `slice_taxa_pq` | `arrange_taxa_pq` | `rename_taxa_pq` |
| occurrences | `filter_occurrences_pq` | - | `mutate_occurrences_pq` | - | - | - |
| tree | `filter_tree_pq` | - | - | - | - | - |

### Data Masking Pattern

Functions use rlang data masking with a special `.` pronoun that refers to the phyloseq object:

```r
filter_samples_pq(data_fungi, Height == "Low", sample_sums(.) > 1000)
```

## Coding Conventions

- Use base pipe (`|>`) not magrittr (`%>%`)
- Use `\() ...` for single-line anonymous functions, `function() {...}` otherwise
- Line length limit: 120 characters
- Tests for `R/{name}.R` go in `tests/testthat/test-{name}.R`

## Documentation Requirements

- Every user-facing function must be exported with roxygen2 documentation
- Wrap roxygen comments at 80 characters
- Parameter docs start with type: `@param x (Int, required) An integer vector of length 1`

## Test Writing Guidelines

- Use `testthat` framework
- Name test files `test-{name}.R` for `R/{name}.R`
- Each test case should have a descriptive name
- Cover edge cases and typical usage scenarios
