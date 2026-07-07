# Changelog

## tidypq 0.2.0 (Development version)

### Breaking changes

- Contamination handling is unified around a detector then filter
  design: each `identify_contam_*_pq()` detector returns a `contam_tbl`
  (one row per flagged taxon), and the single verb
  [`filter_contam_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_contam_pq.md)
  removes the flagged taxa. The previous one-shot functions are removed
  with no deprecation, so old scripts fail with
  `could not find function`. Update calls as follows: `contam_corr_pq()`
  to
  [`identify_contam_corr_pq()`](https://adrientaudiere.github.io/tidypq/reference/identify_contam_corr_pq.md);
  `contam_blocklist_pq()` to
  [`identify_contam_blocklist_pq()`](https://adrientaudiere.github.io/tidypq/reference/identify_contam_blocklist_pq.md);
  `identify_taxa_primer_pq()` and `filter_taxa_primer_pq()` to
  [`identify_contam_primer_pq()`](https://adrientaudiere.github.io/tidypq/reference/identify_contam_primer_pq.md);
  `chimera_removal_dada2()` and `chimera_removal_vs_ref()` to
  `identify_contam_chimera_pq(method = "dada2")` and
  `identify_contam_chimera_pq(method = "vsearch_ref")`;
  `neg_control_classify_pq()`, `neg_control_clean_pq()`, and
  `neg_control_diag_pq()` to
  [`identify_contam_negcontrol_pq()`](https://adrientaudiere.github.io/tidypq/reference/identify_contam_negcontrol_pq.md).
  In every case, pass the result to
  [`filter_contam_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_contam_pq.md)
  to remove taxa.
- The rich diagnostic plots formerly bundled in detector return values
  are removed; a lightweight
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) /
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  method on `contam_tbl` replaces them, with publication-grade figures
  deferred to the `ggplotpq` package.
- `decontam_sam_control()` and `decontam_taxa_control()` are renamed
  [`decontam_control_samples_pq()`](https://adrientaudiere.github.io/tidypq/reference/decontam_control_samples_pq.md)
  and
  [`decontam_control_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/decontam_control_taxa_pq.md)
  (old names fail with `could not find function`); their behaviour is
  unchanged. They remain control-based per-occurrence cleaners and are
  deliberately outside the `contam_tbl` detector family.

### New features

- [`uncross2_pq()`](https://adrientaudiere.github.io/tidypq/reference/uncross2_pq.md)
  can now estimate the rate of tax-switching using internal function
  [`estimate_xtalk_rate_pq()`](https://adrientaudiere.github.io/tidypq/reference/estimate_xtalk_rate_pq.md).
  It estimates the cross-talk (tag-jump) rate `f` de novo from an OTU
  table following Edgar 2018 (UNCROSS2, eqs. 1-2), with
  `rate = "edgar2018"` (paper eq. 2) or `rate = "fraction"` (the raw
  low-read fraction that usearch reports), returning `NA` when the
  estimate is undetermined.
- [`filter_chimera_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_chimera_pq.md)
  detects and removes chimeric taxa in one step (a convenience wrapper
  over
  [`identify_contam_chimera_pq()`](https://adrientaudiere.github.io/tidypq/reference/identify_contam_chimera_pq.md)
  and
  [`filter_contam_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_contam_pq.md)),
  supporting `method = "dada2"` and `method = "vsearch_ref"`.
- [`filter_contam_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_contam_pq.md)
  removes the taxa flagged in a `contam_tbl` (from one detector, or
  several combined with [`rbind()`](https://rdrr.io/r/base/cbind.html)),
  returning the phyloseq object unchanged when nothing is flagged.
- [`merge_markers_pq()`](https://adrientaudiere.github.io/tidypq/reference/merge_markers_pq.md)
  combines several phyloseq objects sequenced with different markers
  (e.g. 16S and ITS) into a single object: it full-joins the `sam_data`,
  disambiguates taxa names by their marker so taxa from different
  objects are never merged, records the marker of origin in a new
  `tax_table` column, and zero-fills unobserved marker/sample
  occurrences.
- [`tax_table_to_df()`](https://adrientaudiere.github.io/tidypq/reference/tax_table_to_df.md)
  converts a phyloseq `tax_table` into a tibble with numeric and list
  columns, lifting the character-only limitation of the `tax_table`
  matrix: numeric-looking ranks are converted to their natural type via
  `convert`, and columns named in `split` become genuine list columns.
- [`identify_contam_blocklist_pq()`](https://adrientaudiere.github.io/tidypq/reference/identify_contam_blocklist_pq.md)
  flags taxa whose genus matches a curated blocklist of known reagent
  and laboratory contaminants (Salter et al. 2014, micRoclean, GRIMER);
  [`known_contaminant_genera()`](https://adrientaudiere.github.io/tidypq/reference/known_contaminant_genera.md)
  returns the default blocklist.
- [`identify_contam_chimera_pq()`](https://adrientaudiere.github.io/tidypq/reference/identify_contam_chimera_pq.md)
  flags chimeric taxa using `dada2` de novo detection
  (`method = "dada2"`) or vsearch reference-based UCHIME
  (`method = "vsearch_ref"`).
- [`identify_contam_corr_pq()`](https://adrientaudiere.github.io/tidypq/reference/identify_contam_corr_pq.md)
  flags taxa whose relative abundance correlates negatively with sample
  depth (GRIMER-inspired), with `cor_abs`, `cor_rel`, `slope`, and
  `p_value` evidence columns.
- [`identify_contam_negcontrol_pq()`](https://adrientaudiere.github.io/tidypq/reference/identify_contam_negcontrol_pq.md)
  flags taxa from negative-control occurrence patterns and records a
  `subtype` of `artifact`, `lab_contaminant`, or `sample_contaminant`;
  only `artifact` and `lab_contaminant` are flagged by default.
- [`identify_contam_primer_pq()`](https://adrientaudiere.github.io/tidypq/reference/identify_contam_primer_pq.md)
  flags taxa whose representative sequence embeds a primer (forward or
  reverse complement) using IUPAC-aware matching via Biostrings.
- [`new_contam_tbl()`](https://adrientaudiere.github.io/tidypq/reference/contam_tbl.md),
  [`is_contam_tbl()`](https://adrientaudiere.github.io/tidypq/reference/contam_tbl.md),
  and [`plot()`](https://rdrr.io/r/graphics/plot.default.html) /
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  methods support the `contam_tbl` class returned by every detector.
- [`uncross2_pq()`](https://adrientaudiere.github.io/tidypq/reference/uncross2_pq.md)
  removes per-cell cross-talk (tag-jumps) using Edgar’s published
  UNCROSS2 score (eq. 3), with a fixed rate `f` or `f = "auto"` for the
  de novo estimate, `method = "cut"` or `"subtract"`, and an optional
  `return_scores` for the score and flag matrices; it is a faithful,
  gentler implementation of the paper and does not reproduce the more
  aggressive usearch binary filtering.
- [`filter_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_taxa_pq.md),
  [`mutate_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/mutate_taxa_pq.md),
  [`slice_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/slice_taxa_pq.md),
  [`arrange_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/arrange_taxa_pq.md),
  [`select_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/select_taxa_pq.md),
  and
  [`rename_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/rename_taxa_pq.md)
  now align `tax_table` to
  [`taxa_names()`](https://rdrr.io/pkg/phyloseq/man/taxa_names-methods.html)
  before any positional operation. Previously, phyloseq objects whose
  `otu_table` and `tax_table` held the same taxa in different orders
  (e.g. produced by
  [`MiscMetabar::mumu_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/mumu_pq.html))
  caused the wrong taxa to be selected/mutated, so
  `filter_taxa_pq(pq, Kingdom == "Archaea")` could return Bacteria.
- Every verb that returns a phyloseq object (`filter_*_pq()`,
  `select_*_pq()`, `mutate_*_pq()`, `slice_*_pq()`, `arrange_*_pq()`,
  `rename_*_pq()`,
  [`filter_occurrences_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_occurrences_pq.md),
  [`filter_contam_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_contam_pq.md),
  [`filter_chimera_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_chimera_pq.md),
  [`decontam_control_samples_pq()`](https://adrientaudiere.github.io/tidypq/reference/decontam_control_samples_pq.md),
  [`decontam_control_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/decontam_control_taxa_pq.md))
  now reorders a desynchronised `sam_data` or `tax_table` back to the
  canonical `otu_table` order (emitting a message when a reorder was
  needed) and verifies the result with
  `MiscMetabar::verify_pq(check_order = TRUE)` before returning. This
  guarantees the returned object keeps
  [`sample_names()`](https://rdrr.io/pkg/phyloseq/man/sample_names-methods.html)/[`taxa_names()`](https://rdrr.io/pkg/phyloseq/man/taxa_names-methods.html)
  aligned across slots, fixing a class of silent positional bugs
  (e.g. `sam_data` columns no longer ending up attached to the wrong
  samples).

## tidypq 0.1.1

- `chimera_removal_vs_ref()` removes chimeric sequences using UCHIME
  reference-based detection via vsearch.
- `contam_blocklist_pq()` flags taxa whose genus belongs to a curated
  blocklist of known reagent and laboratory contaminants (compiled from
  Salter et al. 2014, micRoclean and GRIMER), a taxonomy-only method
  requiring no negative controls; the accompanying
  [`known_contaminant_genera()`](https://adrientaudiere.github.io/tidypq/reference/known_contaminant_genera.md)
  returns the default blocklist for inspection or extension.
- `contam_corr_pq()` detects potential contaminant taxa by examining the
  correlation between each taxon’s relative abundance and total sample
  reads. Contaminants show a negative correlation (higher proportion in
  low-biomass samples). Returns per-taxon statistics (`cor_abs`,
  `cor_rel`, `slope`, `p_value`), a flagged contaminants subset, and a
  three-panel diagnostic plot. Inspired by the GRIMER tool.
- `neg_control_classify_pq()` classifies taxa found in negative controls
  into `artifact`, `lab_contaminant` and `sample_contaminant` categories
  based on read abundance and occurrence patterns, returning the
  per-taxon classification together with summary, distribution, heatmap
  and ordination plots.

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
- `decontam_sam_control()` and `decontam_taxa_control()` remove
  contaminants using negative controls or control taxa.
- `chimera_removal_dada2()` removes chimeric sequences using ‘dada2’.
- [`create_chimera_pq()`](https://adrientaudiere.github.io/tidypq/reference/create_chimera_pq.md)
  generates synthetic chimeric sequences for benchmarking chimera
  detection.
