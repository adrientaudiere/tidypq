# Detector: taxa whose representative sequence embeds a primer

################################################################################
#' Detect taxa whose representative sequence contains a primer
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' For each taxon, checks whether the representative sequence in the `@@refseq`
#' slot matches any primer in `primers` (forward or reverse complement) using
#' IUPAC-aware pattern matching. Taxa embedding a primer are likely chimeric or
#' otherwise artefactual.
#'
#' This is a pure detector: it flags taxa and returns a [contam_tbl] without
#' modifying `physeq`. Use [filter_contam_pq()] to remove the flagged taxa.
#'
#' @param physeq (phyloseq, required) A phyloseq object with a populated
#'   `@@refseq` slot containing DNA sequences.
#' @param primers (named character vector, required) Primer sequences to search
#'   for. IUPAC ambiguity codes are supported (e.g. `M`, `R`, `Y`, `N`). Names
#'   are reported in the `matched_primers` evidence column.
#' @param verbose (logical, default `TRUE`) Print a one-line detection summary.
#'
#' @return A [contam_tbl] with one row per flagged taxon, `method = "primer"`,
#'   and the evidence columns `matched_primers` (comma-separated primer names)
#'   and `n_reads`. Empty when no taxon embeds a primer.
#' @export
#' @author Adrien Taudière
#'
#' @seealso [filter_contam_pq()] to remove flagged taxa,
#'   [MiscMetabar::cutadapt_remove_primers()] for cutadapt-based primer removal
#'   of fastq files.
#'
#' @examplesIf rlang::is_installed("Biostrings")
#' library(MiscMetabar)
#' # Requires a phyloseq object with a populated @@refseq slot.
#' primers <- c(
#'   fwd = "CCCTACGGGGTGCASCAG",
#'   rev = "GGACTACVSGGGTATCTAAT"
#' )
#' identify_contam_primer_pq(data_fungi, primers)
identify_contam_primer_pq <- function(physeq, primers, verbose = TRUE) {
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

  match_mat <- vapply(
    names(primers),
    function(p_name) {
      p <- Biostrings::DNAString(primers[[p_name]])
      rc_p <- Biostrings::reverseComplement(p)
      fwd <- Biostrings::vcountPattern(p, refseq, fixed = FALSE) > 0
      rev <- Biostrings::vcountPattern(rc_p, refseq, fixed = FALSE) > 0
      fwd | rev
    },
    logical(length(refseq))
  )

  any_match <- rowSums(match_mat) > 0
  bad_idx <- which(any_match)

  flagged <- tibble::tibble(
    taxon = names(refseq)[bad_idx],
    method = "primer",
    matched_primers = vapply(
      bad_idx,
      function(i) {
        paste(names(primers)[match_mat[i, ]], collapse = ", ")
      },
      character(1)
    ),
    n_reads = phyloseq::taxa_sums(physeq)[bad_idx]
  )

  if (verbose) {
    cli::cli_alert_info(
      "primer: flagged {nrow(flagged)} of {length(refseq)} taxa embedding a primer."
    )
  }

  new_contam_tbl(flagged)
}
