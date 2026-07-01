# Removal seam for contaminant taxa flagged by identify_contam_*_pq() detectors

################################################################################
#' Remove contaminant taxa flagged by a detector
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' The single removal verb for contamination handling. It consumes any
#' `contam_tbl` produced by an `identify_contam_*_pq()` detector (or several,
#' combined with [rbind()]) and prunes the flagged taxa from `physeq` by name.
#' When the `contam_tbl` is empty, `physeq` is returned unchanged.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param contam (contam_tbl, required) Output of an `identify_contam_*_pq()`
#'   detector. Removal is keyed on the `taxon` column, so taxa are matched by
#'   name rather than position.
#' @param clean_phyloseq_object (logical, default `TRUE`) Whether to clean the
#'   resulting object with `MiscMetabar::clean_pq()` to drop emptied
#'   taxa/samples.
#' @param verbose (logical, default `TRUE`) Whether to report how many taxa
#'   were removed.
#'
#' @return A phyloseq object with the flagged taxa removed.
#' @export
#' @author Adrien Taudière
#'
#' @seealso [identify_contam_corr_pq()] and the other `identify_contam_*_pq()`
#'   detectors that produce a `contam_tbl`.
#'
#' @examples
#' library(MiscMetabar)
#' \donttest{
#' flagged <- identify_contam_corr_pq(data_fungi, verbose = FALSE)
#' data_clean <- filter_contam_pq(data_fungi, flagged)
#' }
filter_contam_pq <- function(
  physeq,
  contam,
  clean_phyloseq_object = TRUE,
  verbose = TRUE
) {
  MiscMetabar::verify_pq(physeq, check_order = FALSE)
  physeq <- canonicalize_pq_order(physeq)
  if (!is_contam_tbl(contam)) {
    stop(
      "`contam` must be a contam_tbl, the output of an ",
      "identify_contam_*_pq() detector.",
      call. = FALSE
    )
  }

  bad <- unique(contam$taxon)
  if (length(bad) == 0) {
    if (verbose) {
      cli::cli_alert_success("No contaminant taxa to remove.")
    }
    return(physeq)
  }

  keep <- !phyloseq::taxa_names(physeq) %in% bad
  if (verbose) {
    cli::cli_alert_warning(
      "Removing {sum(!keep)} contaminant taxon{?s} flagged by method{?s} {.val {unique(contam$method)}}."
    )
  }

  res <- phyloseq::prune_taxa(keep, physeq)
  if (clean_phyloseq_object) {
    res <- MiscMetabar::clean_pq(res, silent = TRUE)
  }
  MiscMetabar::verify_pq(res)
  res
}
