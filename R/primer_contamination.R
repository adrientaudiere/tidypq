# Primer contamination detection and removal

################################################################################
#' Detect taxa whose representative sequence contains a primer
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' For each taxon in `physeq`, checks whether the representative sequence in
#' the `@@refseq` slot matches any primer from `primers` (forward or reverse
#' complement) using IUPAC-aware pattern matching (`fixed = FALSE`). Taxa whose
#' sequences embed a primer are likely chimeric or otherwise artefactual and
#' should be inspected or removed.
#'
#' @param physeq (phyloseq, required) A phyloseq object with a populated
#'   `@@refseq` slot containing DNA sequences.
#' @param primers (named character vector, required) Primer sequences to search
#'   for. IUPAC ambiguity codes are supported (e.g. `M`, `R`, `Y`, `N`).
#'   Names are used in the `matched_primers` output column.
#'
#' @return A `data.frame` with one row per contaminated taxon and three columns:
#'   \describe{
#'     \item{`taxon`}{Taxon name (from `taxa_names(physeq)`).}
#'     \item{`matched_primers`}{Comma-separated names of all primers that
#'       matched the sequence (forward or reverse complement).}
#'     \item{`n_reads`}{Total reads across all samples (`taxa_sums(physeq)`).}
#'   }
#'   Returns `NULL` invisibly when no contaminated taxa are found.
#'
#' @export
#' @author Adrien Taudière
#'
#' @seealso [filter_taxa_primer_pq()] to remove the detected taxa,
#'   [contam_corr_pq()] for correlation-based contaminant detection,
#'   [MiscMetabar::cutadapt_remove_primers()] for cutadapt-based contaminant
#'   removal of fastq files.
#'
#' @examples
#' library(MiscMetabar)
#' # Requires a phyloseq object with a refseq slot and Biostrings installed.
#' primers <- c(
#'   fwd = "CCCTACGGGGTGCASCAG",
#'   rev = "GGACTACVSGGGTATCTAAT"
#' )
#' identify_taxa_primer_pq(data_fungi, primers)
identify_taxa_primer_pq <- function(physeq, primers) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop(
      "Package 'Biostrings' is required. Install with: BiocManager::install('Biostrings')",
      call. = FALSE
    )
  }
  MiscMetabar::verify_pq(physeq)
  if (is.null(physeq@refseq)) {
    stop("physeq must have a populated @refseq slot.", call. = FALSE)
  }
  if (is.null(names(primers)) || any(names(primers) == "")) {
    stop("primers must be a named character vector.", call. = FALSE)
  }

  refseq <- physeq@refseq

  match_mat <- vapply(names(primers), function(p_name) {
    p <- Biostrings::DNAString(primers[[p_name]])
    rc_p <- Biostrings::reverseComplement(p)
    fwd <- Biostrings::vcountPattern(p, refseq, fixed = FALSE) > 0
    rev <- Biostrings::vcountPattern(rc_p, refseq, fixed = FALSE) > 0
    fwd | rev
  }, logical(length(refseq)))

  any_match <- rowSums(match_mat) > 0
  if (!any(any_match)) {
    return(invisible(NULL))
  }

  bad_idx <- which(any_match)
  data.frame(
    taxon = names(refseq)[bad_idx],
    matched_primers = apply(
      match_mat[bad_idx, , drop = FALSE], 1,
      function(row) paste(names(primers)[row], collapse = ", ")
    ),
    n_reads = phyloseq::taxa_sums(physeq)[bad_idx],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

################################################################################
#' Remove taxa whose representative sequence contains a primer
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Calls [identify_taxa_primer_pq()] to identify taxa whose representative
#' sequence in `@@refseq` embeds any primer (forward or reverse complement),
#' then removes them with `phyloseq::prune_taxa()`. If no contaminated taxa
#' are found, the original object is returned unchanged.
#'
#' @param physeq (phyloseq, required) A phyloseq object with a populated
#'   `@@refseq` slot containing DNA sequences.
#' @param primers (named character vector, required) Primer sequences to search
#'   for. IUPAC ambiguity codes are supported (e.g. `M`, `R`, `Y`, `N`).
#'   Names are used in the diagnostic message.
#' @param verbose (logical, default `TRUE`) If `TRUE`, prints a message
#'   reporting the number of taxa and reads removed.
#'
#' @return A phyloseq object with primer-contaminated taxa removed.
#'
#' @export
#' @author Adrien Taudière
#'
#' @seealso [identify_taxa_primer_pq()] to inspect contaminated taxa before
#'   removal, [contam_corr_pq()] for correlation-based contaminant detection,
#'   [MiscMetabar::cutadapt_remove_primers()] for cutadapt-based contaminant
#'   removal of fastq files.
#'
#' @examples
#' library(MiscMetabar)
#' # Requires a phyloseq object with a refseq slot and Biostrings installed.
#' primers <- c(
#'   fwd = "CCCTACGGGGT",
#'   rev = "GGTATCTAAT"
#' )
#' # Inspect first
#' identify_taxa_primer_pq(data_fungi, primers)
#' # Then remove
#' data_f_clean <- filter_taxa_primer_pq(data_fungi, primers)

filter_taxa_primer_pq <- function(physeq, primers, verbose = TRUE) {
  bad <- identify_taxa_primer_pq(physeq, primers)
  if (is.null(bad) || nrow(bad) == 0) {
    if (verbose) {
      cli::cli_alert_success("No taxa with primer contamination found.")
    }
    return(physeq)
  }
  if (verbose) {
    cli::cli_alert_warning(
      "Removing {nrow(bad)} taxa containing primer sequences ({sum(bad$n_reads)} reads)."
    )
  }
  phyloseq::prune_taxa(
    !phyloseq::taxa_names(physeq) %in% bad$taxon,
    physeq
  )
}
