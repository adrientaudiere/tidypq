# tidypq 0.3.0 (Development version)

# tidypq 0.2.0
## Breaking changes

* Contamination handling is unified around a detector then filter design: each `identify_contam_*_pq()` detector returns a `contam_tbl` (one row per flagged taxon), and the single verb `filter_contam_pq()` removes the flagged taxa. The previous one-shot functions are removed with no deprecation, so old scripts fail with `could not find function`. Update calls as follows: `contam_corr_pq()` to `identify_contam_corr_pq()`; `contam_blocklist_pq()` to `identify_contam_blocklist_pq()`; `identify_taxa_primer_pq()` and `filter_taxa_primer_pq()` to `identify_contam_primer_pq()`; `chimera_removal_dada2()` and `chimera_removal_vs_ref()` to `identify_contam_chimera_pq(method = "dada2")` and `identify_contam_chimera_pq(method = "vsearch_ref")`; `neg_control_classify_pq()`, `neg_control_clean_pq()`, and `neg_control_diag_pq()` to `identify_contam_negcontrol_pq()`. In every case, pass the result to `filter_contam_pq()` to remove taxa.
* The rich diagnostic plots formerly bundled in detector return values are removed; a lightweight `plot()` / `autoplot()` method on `contam_tbl` replaces them, with publication-grade figures deferred to the `ggplotpq` package.
* `decontam_sam_control()` and `decontam_taxa_control()` are renamed `decontam_control_samples_pq()` and `decontam_control_taxa_pq()` (old names fail with `could not find function`); their behaviour is unchanged. They remain control-based per-occurrence cleaners and are deliberately outside the `contam_tbl` detector family.

## New features

* `uncross2_pq()` can now estimate the rate of tax-switching using internal function `estimate_xtalk_rate_pq()`. It estimates the cross-talk (tag-jump) rate `f` de novo from an OTU table following Edgar 2018 (UNCROSS2, eqs. 1-2), with `rate = "edgar2018"` (paper eq. 2) or `rate = "fraction"` (the raw low-read fraction that usearch reports), returning `NA` when the estimate is undetermined.
* `filter_chimera_pq()` detects and removes chimeric taxa in one step (a convenience wrapper over `identify_contam_chimera_pq()` and `filter_contam_pq()`), supporting `method = "dada2"` and `method = "vsearch_ref"`.
* `filter_contam_pq()` removes the taxa flagged in a `contam_tbl` (from one detector, or several combined with `rbind()`), returning the phyloseq object unchanged when nothing is flagged.
* `merge_markers_pq()` combines several phyloseq objects sequenced with different markers (e.g. 16S and ITS) into a single object: it full-joins the `sam_data`, disambiguates taxa names by their marker so taxa from different objects are never merged, records the marker of origin in a new `tax_table` column, and zero-fills unobserved marker/sample occurrences.
* `tax_table_to_df()` converts a phyloseq `tax_table` into a tibble with numeric and list columns, lifting the character-only limitation of the `tax_table` matrix: numeric-looking ranks are converted to their natural type via `convert`, and columns named in `split` become genuine list columns.
* `identify_contam_blocklist_pq()` flags taxa whose genus matches a curated blocklist of known reagent and laboratory contaminants (Salter et al. 2014, micRoclean, GRIMER); `known_contaminant_genera()` returns the default blocklist.
* `identify_contam_chimera_pq()` flags chimeric taxa using `dada2` de novo detection (`method = "dada2"`) or vsearch reference-based UCHIME (`method = "vsearch_ref"`).
* `identify_contam_corr_pq()` flags taxa whose relative abundance correlates negatively with sample depth (GRIMER-inspired), with `cor_abs`, `cor_rel`, `slope`, and `p_value` evidence columns.
* `identify_contam_negcontrol_pq()` flags taxa from negative-control occurrence patterns and records a `subtype` of `artifact`, `lab_contaminant`, or `sample_contaminant`; only `artifact` and `lab_contaminant` are flagged by default.
* `identify_contam_primer_pq()` flags taxa whose representative sequence embeds a primer (forward or reverse complement) using IUPAC-aware matching via Biostrings.
* `new_contam_tbl()`, `is_contam_tbl()`, and `plot()` / `autoplot()` methods support the `contam_tbl` class returned by every detector.
* `uncross2_pq()` removes per-cell cross-talk (tag-jumps) using Edgar's published UNCROSS2 score (eq. 3), with a fixed rate `f` or `f = "auto"` for the de novo estimate, `method = "cut"` or `"subtract"`, and an optional `return_scores` for the score and flag matrices; it is a faithful, gentler implementation of the paper and does not reproduce the more aggressive usearch binary filtering.
* `filter_taxa_pq()`, `mutate_taxa_pq()`, `slice_taxa_pq()`, `arrange_taxa_pq()`, `select_taxa_pq()`, and `rename_taxa_pq()` now align `tax_table` to `taxa_names()` before any positional operation. Previously, phyloseq objects whose `otu_table` and `tax_table` held the same taxa in different orders (e.g. produced by `MiscMetabar::mumu_pq()`) caused the wrong taxa to be selected/mutated, so `filter_taxa_pq(pq, Kingdom == "Archaea")` could return Bacteria.
* Every verb that returns a phyloseq object (`filter_*_pq()`, `select_*_pq()`, `mutate_*_pq()`, `slice_*_pq()`, `arrange_*_pq()`, `rename_*_pq()`, `filter_occurrences_pq()`, `filter_contam_pq()`, `filter_chimera_pq()`, `decontam_control_samples_pq()`, `decontam_control_taxa_pq()`) now reorders a desynchronised `sam_data` or `tax_table` back to the canonical `otu_table` order (emitting a message when a reorder was needed) and verifies the result with `MiscMetabar::verify_pq(check_order = TRUE)` before returning. This guarantees the returned object keeps `sample_names()`/`taxa_names()` aligned across slots, fixing a class of silent positional bugs (e.g. `sam_data` columns no longer ending up attached to the wrong samples).

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
