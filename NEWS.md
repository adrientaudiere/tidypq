# tidypq (development version)

* `chimera_removal_vs_ref()` removes chimeric sequences using UCHIME reference-based detection via vsearch.

# tidypq 0.1.0

* Initial CRAN submission.
* `filter_samples_pq()`, `select_samdata_pq()`, `mutate_samdata_pq()`, `slice_samples_pq()`, `arrange_samples_pq()`, and `rename_samples_pq()` provide tidyverse-style verbs for manipulating sample data in 'phyloseq' objects.
* `filter_taxa_pq()`, `select_taxa_pq()`, `mutate_taxa_pq()`, `slice_taxa_pq()`, `arrange_taxa_pq()`, and `rename_taxa_pq()` provide tidyverse-style verbs for manipulating taxonomy data.
* `filter_occurrences_pq()` and `mutate_occurrences_pq()` allow cell-level operations on the OTU table with special variables (`sample_total`, `taxon_total`, etc.).
* `filter_tree_pq()` filters phylogenetic trees by taxa names or tip label patterns.
* `taxa_prevalence()` computes per-taxon prevalence across samples.
* `decontam_sam_control()` and `decontam_taxa_control()` remove contaminants using negative controls or control taxa.
* `chimera_removal_dada2()` removes chimeric sequences using 'dada2'.
* `create_chimera_pq()` generates synthetic chimeric sequences for benchmarking chimera detection.
* `plot_sample_depth_pq()` produces diagnostic plots for detecting sequencing depth outliers.
