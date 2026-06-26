# `decontam_control_*_pq()` are control-based decontamination, not detectors

The two control-based functions are renamed `decontam_control_samples_pq()`
(background estimated from negative/blank control *samples*, per taxon) and
`decontam_control_taxa_pq()` (background estimated from spike-in/control *taxa*,
per sample), and are kept deliberately **outside** the `identify_contam_*_pq()`
→ `filter_contam_pq()` detector family (ADR-0001).

They do not flag whole taxa for removal and so cannot honestly produce a
`contam_tbl`: each estimates a noise threshold from a designated control set and
**zeroes individual occurrences** at or below it (and `decontam_control_taxa_pq`
can additionally drop the control taxa). Folding them into the detector pattern
would change their semantics from per-cell subtraction to whole-taxon removal.

The naming axis is the **control reference** (samples vs taxa), not the data
scale — an earlier "occurrence-scale" framing was rejected as misleading,
because the filter metric is a control-derived threshold rather than the
occurrence value itself. `known_contaminant_genera()` and `create_chimera_pq()`
are likewise retained outside the detector family (a blocklist data source and a
synthetic-chimera generator, respectively).
