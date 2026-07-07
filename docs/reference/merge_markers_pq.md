# Merge phyloseq objects from different markers into a single object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Combine several phyloseq objects that describe the **same set of
samples** sequenced with **different markers** (e.g. 16S and ITS, or
ITS1 and ITS2) into a single phyloseq object. The result is assembled
slot by slot with full control (it does *not* rely on
[`phyloseq::merge_phyloseq()`](https://rdrr.io/pkg/phyloseq/man/merge_phyloseq.html),
which merges phyloseq *components* rather than whole objects):

- **Samples** are combined with a *full join* of the `sam_data` slots:
  the output keeps the union of all sample names and the union of all
  metadata columns. A sample shared by two markers keeps a single row;
  for a shared column, the first non-missing value across the input
  objects is used.

- **Taxa are never merged.** Taxa names are disambiguated by prefixing
  them with the marker name (see `disambiguate_taxa`), so a taxon named
  `"ASV_1"` in two markers becomes e.g. `"16S_ASV_1"` and `"ITS_ASV_1"`
  and stays as two distinct taxa. The output `otu_table` is therefore
  the block-diagonal union of the inputs.

- A **new `tax_table` column** (named by `marker_column`) records the
  marker of origin for every taxon.

- An **occurrence** of a taxon in a sample that was not sequenced with
  that taxon's marker is set to `0` (the taxon was simply not observable
  there).

- `refseq` slots are concatenated when *every* input object has one;
  otherwise `refseq` is dropped with a warning. Phylogenetic trees
  cannot be meaningfully merged across markers and are dropped with a
  warning.

## Usage

``` r
merge_markers_pq(
  listpq,
  marker_names = names(listpq),
  marker_column = "marker",
  sep = "_",
  disambiguate_taxa = TRUE
)
```

## Arguments

- listpq:

  (named list of phyloseq, required) A list of at least two phyloseq
  objects, one per marker. If the list is named, the names are used as
  marker names (see `marker_names`).

- marker_names:

  (character, default `names(listpq)`) Marker name for each object, in
  the same order as `listpq`. Must be unique and have one entry per
  object. Required when `listpq` has no names.

- marker_column:

  (character of length 1, default `"marker"`) Name of the new
  `tax_table` column that stores the marker of origin.

- sep:

  (character of length 1, default `"_"`) Separator inserted between the
  marker name and the original taxon name when disambiguating.

- disambiguate_taxa:

  (logical, default `TRUE`) If `TRUE`, every taxon name is prefixed with
  its marker name to guarantee that no two taxa from different objects
  can ever be merged. If `FALSE`, the original taxa names are kept and
  the function errors when names collide across objects.

## Value

A single phyloseq object stacking all input markers (union of samples,
disambiguated union of taxa, marker column in `tax_table`).

## Do NOT post-cluster the result

The taxa carried by different markers are intentionally kept apart. You
must **not** run any post-clustering step (e.g.
[`MiscMetabar::postcluster_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.html))
on the merged object, because that would collapse taxa coming from
different markers into a single taxon, mixing incomparable sequences (a
16S OTU and an ITS OTU are not the same organism just because their
sequences happen to be similar). Keep the markers separate for any
per-marker analysis and treat the merged object as a *stacked* view for
cross-marker exploration.

The only legitimate cross-marker merging is at the **taxonomy** level
and must be done explicitly on the `tax_table` (e.g. aggregating ITS1
and ITS2 taxa that resolve to the same species with
[`MiscMetabar::merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.html)
or
[`phyloseq::tax_glom()`](https://rdrr.io/pkg/phyloseq/man/tax_glom.html)),
never through sequence-based clustering.

## Author

Adrien Taudière

## Examples

``` r
library(phyloseq)

# Two tiny fake markers describing the same three samples ------------------
make_marker <- function(taxa, ranks) {
  otu <- matrix(
    sample(0:20, length(taxa) * 3, replace = TRUE),
    nrow = length(taxa),
    dimnames = list(taxa, c("samp1", "samp2", "samp3"))
  )
  tax <- matrix(
    paste(rep(ranks, each = length(taxa)), seq_along(taxa), sep = "_"),
    nrow = length(taxa),
    dimnames = list(taxa, ranks)
  )
  sam <- data.frame(
    row.names = c("samp1", "samp2", "samp3"),
    site = c("A", "A", "B"),
    depth = c(10, 20, 30)
  )
  phyloseq(
    otu_table(otu, taxa_are_rows = TRUE),
    tax_table(tax),
    sample_data(sam)
  )
}

pq_16S <- make_marker(c("ASV_1", "ASV_2"), c("Kingdom", "Phylum"))
# Note the SAME taxon name "ASV_1" reused in the ITS object:
pq_ITS <- make_marker(c("ASV_1", "ASV_3"), c("Kingdom", "Phylum", "Genus"))

merged <- merge_markers_pq(list("marker_16S" = pq_16S, "marker_ITS" = pq_ITS))
#> Merged 2 markers (marker_16S, marker_ITS) into 4 taxa across 3 samples.

# Taxa are disambiguated and never merged: 4 taxa, not 3.
taxa_names(merged)
#> [1] "marker_16S_ASV_1" "marker_16S_ASV_2" "marker_ITS_ASV_1" "marker_ITS_ASV_3"
# The marker of origin is recorded in a new tax_table column.
tax_table(merged)[, "marker"]
#> Taxonomy Table:     [4 taxa by 1 taxonomic ranks]:
#>                  marker      
#> marker_16S_ASV_1 "marker_16S"
#> marker_16S_ASV_2 "marker_16S"
#> marker_ITS_ASV_1 "marker_ITS"
#> marker_ITS_ASV_3 "marker_ITS"

# Providing marker names explicitly for an unnamed list:
merge_markers_pq(list(pq_16S, pq_ITS), marker_names = c("16S", "ITS"))
#> Merged 2 markers (16S, ITS) into 4 taxa across 3 samples.
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 4 taxa and 3 samples ]
#> sample_data() Sample Data:       [ 3 samples by 2 sample variables ]
#> tax_table()   Taxonomy Table:    [ 4 taxa by 4 taxonomic ranks ]

# WARNING: never post-cluster `merged` (e.g. MiscMetabar::postcluster_pq()).
# That would wrongly merge a 16S taxon with an ITS taxon.
```
