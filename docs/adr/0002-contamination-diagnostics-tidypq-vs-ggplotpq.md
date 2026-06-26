# Contamination diagnostics: lightweight in tidypq, rich figures in ggplotpq

Contamination diagnostics are split out of detectors (per ADR-0001) and live in
two places. tidypq ships a lightweight `autoplot()`/`plot()` method on
`contam_tbl` — faceted by `method`, with evidence columns driving the geom —
using the `ggplot2`/`patchwork` it already imports. Publication-grade
contamination figures are deferred to ggplotpq (consuming the same `contam_tbl`)
and built only if demand appears.

This is a deliberate, partial deviation from the pqverse rule that **ggplotpq
owns visualisation**. A reader who finds plotting code in tidypq should know it
is intentional: a contamination diagnostic is tightly coupled to its detector's
evidence columns (a correlation scatter needs `cor`/`slope`; a blocklist match
has nothing to plot), so keeping the basic view next to the detector buys
locality and adds no new dependency. The `contam_tbl` contract is the seam that
lets richer figures move to ggplotpq later without changing tidypq.

Realised: the four-panel negative-control diagnostic `neg_control_diag_pq()` was
migrated out of tidypq into `ggplotpq` (its programmatic counterpart stays here
as `identify_contam_negcontrol_pq()`). It is the first rich contamination figure
to live in ggplotpq under this split.
