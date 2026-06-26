# Contamination handling: detector → filter seam with a `contam_tbl` contract

Every contamination method in tidypq is split into a pure **detector**
(`identify_contam_<method>_pq()`) that flags taxa and returns a class-tagged
`contam_tbl` (one row per flagged taxon, keyed by `taxon`, plus a `method`
column and method-specific evidence columns), and a single **act** verb
(`filter_contam_pq()`) that consumes any `contam_tbl` and removes the flagged
taxa. We chose this over the previous pattern, where each function detected,
removed, and plotted in one call and returned an idiosyncratic shape (phyloseq,
data.frame, or a named list bundling stats and plots).

## Considered options

- **Detector → filter seam (chosen).** Detectors compose: `rbind()` two
  `contam_tbl`s and `filter_contam_pq()` removes the union. The interface is
  small (one return type, one removal verb) and deep (each detector hides
  substantial method-specific logic). The `taxon` key avoids the positional
  `tax_table`/`otu_table` misalignment fixed in 0.2.0.
- **Keep contamination as a separate subsystem**, merely made internally
  consistent — rejected: misses the fact that removal *is* taxa-scale filtering.
- **Fat all-in-one functions with a uniform list return** — rejected: optimises
  for one-call convenience at the cost of composability and testability (you
  cannot test detection independently of removal or plotting).

## Consequences

- Released names (`chimera_removal_dada2()`, `contam_corr_pq()`,
  `neg_control_classify_pq()`, `decontam_*_control()`, …) are removed, not
  deprecated: a **hard break** with **bare removal** in the next 0.x, documented
  in `NEWS.md`. Pre-1.0 status makes this acceptable; old scripts fail with R's
  standard "could not find function" rather than a migration message.
- Detection now runs once and yields a reusable `contam_tbl`; expensive
  detectors (chimera via dada2/vsearch) are no longer re-run to re-filter.
- One deliberate exception to the "no one-shot wrapper" stance: `filter_chimera_pq()`
  composes `identify_contam_chimera_pq()` + `filter_contam_pq()` in a single call,
  because chimera removal is a near-ubiquitous cleaning step. No equivalent wrapper
  is provided for the other detectors; they are composed explicitly.
