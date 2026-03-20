# Changelog

## tidypq 0.1.0

- Initial CRAN submission.
- [`filter_samples_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_samples_pq.md),
  [`select_samdata_pq()`](https://adrientaudiere.github.io/tidypq/reference/select_samdata_pq.md),
  [`mutate_samdata_pq()`](https://adrientaudiere.github.io/tidypq/reference/mutate_samdata_pq.md),
  [`slice_samples_pq()`](https://adrientaudiere.github.io/tidypq/reference/slice_samples_pq.md),
  [`arrange_samples_pq()`](https://adrientaudiere.github.io/tidypq/reference/arrange_samples_pq.md),
  and
  [`rename_samples_pq()`](https://adrientaudiere.github.io/tidypq/reference/rename_samples_pq.md)
  provide tidyverse-style verbs for manipulating sample data in
  ‘phyloseq’ objects.
- [`filter_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_taxa_pq.md),
  [`select_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/select_taxa_pq.md),
  [`mutate_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/mutate_taxa_pq.md),
  [`slice_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/slice_taxa_pq.md),
  [`arrange_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/arrange_taxa_pq.md),
  and
  [`rename_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/rename_taxa_pq.md)
  provide tidyverse-style verbs for manipulating taxonomy data.
- [`filter_occurrences_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_occurrences_pq.md)
  and
  [`mutate_occurrences_pq()`](https://adrientaudiere.github.io/tidypq/reference/mutate_occurrences_pq.md)
  allow cell-level operations on the OTU table with special variables
  (`sample_total`, `taxon_total`, etc.).
- [`filter_tree_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_tree_pq.md)
  filters phylogenetic trees by taxa names or tip label patterns.
- [`taxa_prevalence()`](https://adrientaudiere.github.io/tidypq/reference/taxa_prevalence.md)
  computes per-taxon prevalence across samples.
- [`decontam_sam_control()`](https://adrientaudiere.github.io/tidypq/reference/decontam_sam_control.md)
  and
  [`decontam_taxa_control()`](https://adrientaudiere.github.io/tidypq/reference/decontam_taxa_control.md)
  remove contaminants using negative controls or control taxa.
- [`chimera_removal_dada2()`](https://adrientaudiere.github.io/tidypq/reference/chimera_removal_dada2.md)
  removes chimeric sequences using ‘dada2’.
- [`create_chimera_pq()`](https://adrientaudiere.github.io/tidypq/reference/create_chimera_pq.md)
  generates synthetic chimeric sequences for benchmarking chimera
  detection.
- [`plot_sample_depth_pq()`](https://adrientaudiere.github.io/tidypq/reference/plot_sample_depth_pq.md)
  produces diagnostic plots for detecting sequencing depth outliers.
