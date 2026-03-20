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

## Helpers

Utility functions for common operations.

- [`taxa_prevalence()`](https://adrientaudiere.github.io/tidypq/reference/taxa_prevalence.md)
  : Calculate taxa prevalence

## Decontamination

Remove contaminants using negative controls or control taxa.

- [`decontam_sam_control()`](https://adrientaudiere.github.io/tidypq/reference/decontam_sam_control.md)
  : Decontaminate based on negative/blank control samples
- [`decontam_taxa_control()`](https://adrientaudiere.github.io/tidypq/reference/decontam_taxa_control.md)
  : Decontaminate based on control taxa

## Chimera detection

Detect and benchmark chimeric sequences.

- [`chimera_removal_dada2()`](https://adrientaudiere.github.io/tidypq/reference/chimera_removal_dada2.md)
  : Remove chimeric sequences using dada2
- [`create_chimera_pq()`](https://adrientaudiere.github.io/tidypq/reference/create_chimera_pq.md)
  : Create a phyloseq object with synthetic chimeric sequences

## Plots

Diagnostic visualizations.

- [`plot_sample_depth_pq()`](https://adrientaudiere.github.io/tidypq/reference/plot_sample_depth_pq.md)
  : Plot sample depth differences to detect outliers
