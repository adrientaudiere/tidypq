# Merge phyloseq objects issued from different markers

################################################################################
#' Merge phyloseq objects from different markers into a single object
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Combine several phyloseq objects that describe the **same set of samples**
#' sequenced with **different markers** (e.g. 16S and ITS, or ITS1 and ITS2)
#' into a single phyloseq object. The result is assembled slot by slot with
#' full control (it does *not* rely on `phyloseq::merge_phyloseq()`, which
#' merges phyloseq *components* rather than whole objects):
#'
#' - **Samples** are combined with a *full join* of the `sam_data` slots: the
#'   output keeps the union of all sample names and the union of all metadata
#'   columns. A sample shared by two markers keeps a single row; for a shared
#'   column, the first non-missing value across the input objects is used.
#' - **Taxa are never merged.** Taxa names are disambiguated by prefixing them
#'   with the marker name (see `disambiguate_taxa`), so a taxon named `"ASV_1"`
#'   in two markers becomes e.g. `"16S_ASV_1"` and `"ITS_ASV_1"` and stays as
#'   two distinct taxa. The output `otu_table` is therefore the block-diagonal
#'   union of the inputs.
#' - A **new `tax_table` column** (named by `marker_column`) records the marker
#'   of origin for every taxon.
#' - An **occurrence** of a taxon in a sample that was not sequenced with that
#'   taxon's marker is set to `0` (the taxon was simply not observable there).
#' - `refseq` slots are concatenated when *every* input object has one;
#'   otherwise `refseq` is dropped with a warning. Phylogenetic trees cannot be
#'   meaningfully merged across markers and are dropped with a warning.
#'
#' @section Do NOT post-cluster the result:
#' The taxa carried by different markers are intentionally kept apart. You must
#' **not** run any post-clustering step (e.g. `MiscMetabar::postcluster_pq()`)
#' on the merged object, because that would collapse taxa coming from different
#' markers into a single taxon, mixing incomparable sequences (a 16S OTU and an
#' ITS OTU are not the same organism just because their sequences happen to be
#' similar). Keep the markers separate for any per-marker analysis and treat the
#' merged object as a *stacked* view for cross-marker exploration.
#'
#' The only legitimate cross-marker merging is at the **taxonomy** level and
#' must be done explicitly on the `tax_table` (e.g. aggregating ITS1 and ITS2
#' taxa that resolve to the same species with `MiscMetabar::merge_taxa_vec()` or
#' `phyloseq::tax_glom()`), never through sequence-based clustering.
#'
#' @param listpq (named list of phyloseq, required) A list of at least two
#'   phyloseq objects, one per marker. If the list is named, the names are used
#'   as marker names (see `marker_names`).
#' @param marker_names (character, default `names(listpq)`) Marker name for each
#'   object, in the same order as `listpq`. Must be unique and have one entry
#'   per object. Required when `listpq` has no names.
#' @param marker_column (character of length 1, default `"marker"`) Name of the
#'   new `tax_table` column that stores the marker of origin.
#' @param sep (character of length 1, default `"_"`) Separator inserted between
#'   the marker name and the original taxon name when disambiguating.
#' @param disambiguate_taxa (logical, default `TRUE`) If `TRUE`, every taxon
#'   name is prefixed with its marker name to guarantee that no two taxa from
#'   different objects can ever be merged. If `FALSE`, the original taxa names
#'   are kept and the function errors when names collide across objects.
#'
#' @return A single phyloseq object stacking all input markers (union of
#'   samples, disambiguated union of taxa, marker column in `tax_table`).
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' library(phyloseq)
#'
#' # Two tiny fake markers describing the same three samples ------------------
#' make_marker <- function(taxa, ranks) {
#'   otu <- matrix(
#'     sample(0:20, length(taxa) * 3, replace = TRUE),
#'     nrow = length(taxa),
#'     dimnames = list(taxa, c("samp1", "samp2", "samp3"))
#'   )
#'   tax <- matrix(
#'     paste(rep(ranks, each = length(taxa)), seq_along(taxa), sep = "_"),
#'     nrow = length(taxa),
#'     dimnames = list(taxa, ranks)
#'   )
#'   sam <- data.frame(
#'     row.names = c("samp1", "samp2", "samp3"),
#'     site = c("A", "A", "B"),
#'     depth = c(10, 20, 30)
#'   )
#'   phyloseq(
#'     otu_table(otu, taxa_are_rows = TRUE),
#'     tax_table(tax),
#'     sample_data(sam)
#'   )
#' }
#'
#' pq_16S <- make_marker(c("ASV_1", "ASV_2"), c("Kingdom", "Phylum"))
#' # Note the SAME taxon name "ASV_1" reused in the ITS object:
#' pq_ITS <- make_marker(c("ASV_1", "ASV_3"), c("Kingdom", "Phylum", "Genus"))
#'
#' merged <- merge_markers_pq(list("marker_16S" = pq_16S, "marker_ITS" = pq_ITS))
#'
#' # Taxa are disambiguated and never merged: 4 taxa, not 3.
#' taxa_names(merged)
#' # The marker of origin is recorded in a new tax_table column.
#' tax_table(merged)[, "marker"]
#'
#' # Providing marker names explicitly for an unnamed list:
#' merge_markers_pq(list(pq_16S, pq_ITS), marker_names = c("16S", "ITS"))
#'
#' # WARNING: never post-cluster `merged` (e.g. MiscMetabar::postcluster_pq()).
#' # That would wrongly merge a 16S taxon with an ITS taxon.
merge_markers_pq <- function(
  listpq,
  marker_names = names(listpq),
  marker_column = "marker",
  sep = "_",
  disambiguate_taxa = TRUE
) {
  if (!is.list(listpq) || length(listpq) < 2) {
    stop("`listpq` must be a list of at least two phyloseq objects.")
  }
  if (!all(vapply(listpq, \(x) inherits(x, "phyloseq"), logical(1)))) {
    stop("All elements of `listpq` must be phyloseq objects.")
  }
  if (
    is.null(marker_names) ||
      length(marker_names) != length(listpq) ||
      anyNA(marker_names) ||
      any(!nzchar(marker_names))
  ) {
    stop(
      "Provide one non-empty marker name per object, either as the names of ",
      "`listpq` or via the `marker_names` argument."
    )
  }
  if (anyDuplicated(marker_names)) {
    stop("`marker_names` must be unique.")
  }
  if (!is.character(marker_column) || length(marker_column) != 1) {
    stop("`marker_column` must be a single character string.")
  }
  if (!is.character(sep) || length(sep) != 1) {
    stop("`sep` must be a single character string.")
  }

  # Prepare each object: canonical order, disambiguated taxa names, marker
  # column added to the tax_table, taxa-as-rows otu matrix.
  prep <- lapply(seq_along(listpq), function(i) {
    physeq <- listpq[[i]]
    MiscMetabar::verify_pq(physeq, check_order = FALSE)
    physeq <- canonicalize_pq_order(physeq)
    mk <- marker_names[[i]]

    old_taxa <- phyloseq::taxa_names(physeq)
    new_taxa <- if (disambiguate_taxa) {
      paste0(mk, sep, old_taxa)
    } else {
      old_taxa
    }

    otu <- as(MiscMetabar::taxa_as_rows(physeq)@otu_table, "matrix")
    otu <- otu[old_taxa, , drop = FALSE]
    rownames(otu) <- new_taxa

    tax <- taxa_table_df(physeq)
    rownames(tax) <- new_taxa
    if (marker_column %in% colnames(tax)) {
      stop(
        "`marker_column` = '",
        marker_column,
        "' already exists in the tax_table of object ",
        i,
        ". Choose another `marker_column`."
      )
    }
    tax[[marker_column]] <- mk

    sam <- if (!is.null(physeq@sam_data)) {
      as(phyloseq::sample_data(physeq), "data.frame")
    } else {
      NULL
    }

    refseq <- if (!is.null(physeq@refseq)) {
      rs <- physeq@refseq
      rs <- rs[old_taxa]
      names(rs) <- new_taxa
      rs
    } else {
      NULL
    }

    list(
      otu = otu,
      tax = tax,
      sam = sam,
      refseq = refseq,
      has_tree = !is.null(phyloseq::phy_tree(physeq, errorIfNULL = FALSE))
    )
  })

  all_taxa <- unlist(lapply(prep, \(p) rownames(p$otu)), use.names = FALSE)
  if (anyDuplicated(all_taxa)) {
    stop(
      "Taxa names collide across objects. Set `disambiguate_taxa = TRUE` ",
      "(the default) to prefix taxa names with their marker name."
    )
  }
  all_samples <- Reduce(union, lapply(prep, \(p) colnames(p$otu)))

  # otu_table: block union, zero-filled for unobserved (marker, sample) pairs.
  otu_merged <- matrix(
    0,
    nrow = length(all_taxa),
    ncol = length(all_samples),
    dimnames = list(all_taxa, all_samples)
  )
  for (p in prep) {
    otu_merged[rownames(p$otu), colnames(p$otu)] <- p$otu
  }

  # tax_table: union of ranks (NA-filled), marker column pushed last.
  all_tax_cols <- Reduce(union, lapply(prep, \(p) colnames(p$tax)))
  all_tax_cols <- c(setdiff(all_tax_cols, marker_column), marker_column)
  tax_merged <- do.call(
    rbind,
    lapply(prep, function(p) {
      for (col in setdiff(all_tax_cols, colnames(p$tax))) {
        p$tax[[col]] <- NA_character_
      }
      p$tax[, all_tax_cols, drop = FALSE]
    })
  )

  components <- list(
    phyloseq::otu_table(otu_merged, taxa_are_rows = TRUE),
    phyloseq::tax_table(as.matrix(tax_merged))
  )

  # sam_data: full join across objects, reindexed to every merged sample.
  sam_present <- Filter(Negate(is.null), lapply(prep, \(p) p$sam))
  if (length(sam_present) > 0) {
    sam_list <- lapply(sam_present, function(df) {
      df[[".merge_sample"]] <- rownames(df)
      df
    })
    sam_all <- dplyr::bind_rows(sam_list)
    sam_merged <- dplyr::summarise(
      dplyr::group_by(sam_all, .data[[".merge_sample"]]),
      dplyr::across(
        dplyr::everything(),
        \(x) {
          non_na <- x[!is.na(x)]
          if (length(non_na) > 0) non_na[[1]] else x[[1]]
        }
      ),
      .groups = "drop"
    )
    sam_df <- as.data.frame(sam_merged)
    rownames(sam_df) <- sam_df[[".merge_sample"]]
    sam_df[[".merge_sample"]] <- NULL
    sam_df <- sam_df[all_samples, , drop = FALSE]
    rownames(sam_df) <- all_samples
    components <- c(
      components,
      list(phyloseq::sample_data(sam_df))
    )
  }

  has_refseq <- vapply(prep, \(p) !is.null(p$refseq), logical(1))
  if (all(has_refseq)) {
    components <- c(
      components,
      list(do.call(c, lapply(prep, \(p) p$refseq)))
    )
  } else if (any(has_refseq)) {
    warning(
      "`refseq` dropped: it is only present in some of the input objects."
    )
  }

  if (any(vapply(prep, \(p) p$has_tree, logical(1)))) {
    warning(
      "Phylogenetic tree(s) dropped: trees cannot be merged across markers."
    )
  }

  new_physeq <- do.call(phyloseq::phyloseq, components)

  message(
    "Merged ",
    length(listpq),
    " markers (",
    paste(marker_names, collapse = ", "),
    ") into ",
    phyloseq::ntaxa(new_physeq),
    " taxa across ",
    phyloseq::nsamples(new_physeq),
    " samples."
  )

  MiscMetabar::verify_pq(new_physeq)
  return(new_physeq)
}
