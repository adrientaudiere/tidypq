# tidypq 0.2.0 (Development version)

* `filter_taxa_pq()`, `mutate_taxa_pq()`, `slice_taxa_pq()`, `arrange_taxa_pq()`, `select_taxa_pq()`, and `rename_taxa_pq()` now align `tax_table` to `taxa_names()` before any positional operation. Previously, phyloseq objects whose `otu_table` and `tax_table` held the same taxa in different orders (e.g. produced by `MiscMetabar::mumu_pq()`) caused the wrong taxa to be selected/mutated, so `filter_taxa_pq(pq, Kingdom == "Archaea")` could return Bacteria.
* `identify_taxa_primer_pq()` identifies taxa whose representative sequence (`@refseq`) contains any primer (forward or reverse complement) using IUPAC-aware pattern matching via Biostrings. Returns a data.frame of contaminated taxa with matched primer names and read counts.
* `filter_taxa_primer_pq()` removes primer-contaminated taxa detected by `identify_taxa_primer_pq()` and returns the cleaned phyloseq object unchanged when no contamination is found.

# tidypq 0.1.1
* `chimera_removal_vs_ref()` removes chimeric sequences using UCHIME reference-based detection via vsearch.
* `contam_blocklist_pq()` flags taxa whose genus belongs to a curated blocklist of known reagent and laboratory contaminants (compiled from Salter et al. 2014, micRoclean and GRIMER), a taxonomy-only method requiring no negative controls; the accompanying `known_contaminant_genera()` returns the default blocklist for inspection or extension.
* `contam_corr_pq()` detects potential contaminant taxa by examining the correlation between each taxon's relative abundance and total sample reads. Contaminants show a negative correlation (higher proportion in low-biomass samples). Returns per-taxon statistics (`cor_abs`, `cor_rel`, `slope`, `p_value`), a flagged contaminants subset, and a three-panel diagnostic plot. Inspired by the GRIMER tool.
* `neg_control_classify_pq()` classifies taxa found in negative controls into `artifact`, `lab_contaminant` and `sample_contaminant` categories based on read abundance and occurrence patterns, returning the per-taxon classification together with summary, distribution, heatmap and ordination plots.

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
