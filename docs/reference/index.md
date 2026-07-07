# Package index

## Sample-scale verbs

Filter, select, mutate, slice, arrange, and rename samples.

- [`filter_samples_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_samples_pq.md)
  : Filter samples in a phyloseq object
- [`select_samdata_pq()`](https://adrientaudiere.github.io/tidypq/reference/select_samdata_pq.md)
  : Select columns from sample_data in a phyloseq object
- [`mutate_samdata_pq()`](https://adrientaudiere.github.io/tidypq/reference/mutate_samdata_pq.md)
  : Add or modify columns in sample_data
- [`slice_samples_pq()`](https://adrientaudiere.github.io/tidypq/reference/slice_samples_pq.md)
  : Subset samples by position
- [`arrange_samples_pq()`](https://adrientaudiere.github.io/tidypq/reference/arrange_samples_pq.md)
  : Arrange samples by column values
- [`rename_samples_pq()`](https://adrientaudiere.github.io/tidypq/reference/rename_samples_pq.md)
  : Rename columns in sample_data

## Taxa-scale verbs

Filter, select, mutate, slice, arrange, and rename taxa.

- [`filter_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_taxa_pq.md)
  : Filter taxa in a phyloseq object
- [`select_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/select_taxa_pq.md)
  : Select columns from tax_table in a phyloseq object
- [`mutate_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/mutate_taxa_pq.md)
  : Add or modify columns in tax_table
- [`slice_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/slice_taxa_pq.md)
  : Subset taxa by position
- [`arrange_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/arrange_taxa_pq.md)
  : Arrange taxa by column values
- [`rename_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/rename_taxa_pq.md)
  : Rename columns in tax_table

## Occurrence-scale verbs

Filter and transform OTU table values.

- [`filter_occurrences_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_occurrences_pq.md)
  : Filter occurrences in the OTU table
- [`mutate_occurrences_pq()`](https://adrientaudiere.github.io/tidypq/reference/mutate_occurrences_pq.md)
  : Transform OTU table values

## Tree-scale verbs

Filter phylogenetic trees.

- [`filter_tree_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_tree_pq.md)
  : Filter phyloseq by tree topology

## Tidy conversion

Convert a phyloseq object into a tidy tibble.

- [`pq_to_tidy()`](https://adrientaudiere.github.io/tidypq/reference/pq_to_tidy.md)
  : Convert a phyloseq object to a tidy tibble

- [`tax_table_to_df()`](https://adrientaudiere.github.io/tidypq/reference/tax_table_to_df.md)
  :

  Convert a `tax_table` into a data.frame with numeric and list columns

## Merging markers

Stack phyloseq objects from different markers into one object.

- [`merge_markers_pq()`](https://adrientaudiere.github.io/tidypq/reference/merge_markers_pq.md)
  : Merge phyloseq objects from different markers into a single object

## Helpers

Utility functions for common operations.

- [`taxa_prevalence()`](https://adrientaudiere.github.io/tidypq/reference/taxa_prevalence.md)
  : Calculate taxa prevalence

## Contaminant detectors

Flag suspect taxa; each returns a contam_tbl consumed by
filter_contam_pq().

- [`identify_contam_blocklist_pq()`](https://adrientaudiere.github.io/tidypq/reference/identify_contam_blocklist_pq.md)
  : Detect taxa whose genus belongs to a known-contaminant blocklist
- [`identify_contam_chimera_pq()`](https://adrientaudiere.github.io/tidypq/reference/identify_contam_chimera_pq.md)
  : Detect chimeric taxa with dada2 or vsearch
- [`identify_contam_corr_pq()`](https://adrientaudiere.github.io/tidypq/reference/identify_contam_corr_pq.md)
  : Detect contaminants by correlating taxon reads with total sample
  depth
- [`identify_contam_negcontrol_pq()`](https://adrientaudiere.github.io/tidypq/reference/identify_contam_negcontrol_pq.md)
  : Detect contaminants from negative-control occurrence patterns
- [`identify_contam_primer_pq()`](https://adrientaudiere.github.io/tidypq/reference/identify_contam_primer_pq.md)
  : Detect taxa whose representative sequence contains a primer
- [`known_contaminant_genera()`](https://adrientaudiere.github.io/tidypq/reference/known_contaminant_genera.md)
  : Known reagent- and laboratory-contaminant genera

## Contaminant removal and the contam_tbl object

Remove flagged taxa and work with the detector output object.

- [`filter_contam_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_contam_pq.md)
  : Remove contaminant taxa flagged by a detector

- [`filter_chimera_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_chimera_pq.md)
  : Remove chimeric taxa in one step

- [`new_contam_tbl()`](https://adrientaudiere.github.io/tidypq/reference/contam_tbl.md)
  [`is_contam_tbl()`](https://adrientaudiere.github.io/tidypq/reference/contam_tbl.md)
  :

  The `contam_tbl` class

- [`print(`*`<contam_tbl>`*`)`](https://adrientaudiere.github.io/tidypq/reference/print.contam_tbl.md)
  :

  Print a `contam_tbl`

- [`rbind(`*`<contam_tbl>`*`)`](https://adrientaudiere.github.io/tidypq/reference/rbind.contam_tbl.md)
  :

  Combine several `contam_tbl` objects

- [`autoplot(`*`<contam_tbl>`*`)`](https://adrientaudiere.github.io/tidypq/reference/autoplot.contam_tbl.md)
  [`plot(`*`<contam_tbl>`*`)`](https://adrientaudiere.github.io/tidypq/reference/autoplot.contam_tbl.md)
  :

  Diagnostic plot for a `contam_tbl`

## Control-based decontamination

Zero out occurrences below a background level estimated from control
samples or control taxa.

- [`decontam_control_samples_pq()`](https://adrientaudiere.github.io/tidypq/reference/decontam_control_samples_pq.md)
  : Decontaminate based on negative/blank control samples
- [`decontam_control_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/decontam_control_taxa_pq.md)
  : Decontaminate based on control taxa

## Cross-talk (tag-jump) filtering

Detect, visualise, and remove per-cell cross-talk (tag-switching) with
the UNCROSS2 score.

- [`uncross2_pq()`](https://adrientaudiere.github.io/tidypq/reference/uncross2_pq.md)
  : Remove cross-talk (tag-jumps) using the UNCROSS2 score
- [`plot_uncross2_scores()`](https://adrientaudiere.github.io/tidypq/reference/plot_uncross2_scores.md)
  : Plot UNCROSS2 cross-talk scores against OTU abundance

## Chimera utilities

Generate synthetic chimeras for benchmarking detectors.

- [`create_chimera_pq()`](https://adrientaudiere.github.io/tidypq/reference/create_chimera_pq.md)
  : Create a phyloseq object with synthetic chimeric sequences
