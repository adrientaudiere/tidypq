# Detector: chimeric sequences via dada2 (de novo) or vsearch (reference-based)

################################################################################
#' Detect chimeric taxa with dada2 or vsearch
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Flags chimeric taxa using either `dada2::removeBimeraDenovo()` de novo
#' detection (`method = "dada2"`) or UCHIME reference-based detection via
#' [vsearch](https://github.com/torognes/vsearch) (`method = "vsearch_ref"`).
#'
#' This is a pure detector: it flags taxa and returns a [contam_tbl] without
#' modifying `physeq`. Use [filter_contam_pq()] to remove the flagged taxa.
#'
#' @param physeq (phyloseq, required) A phyloseq object with a populated
#'   `@@refseq` slot containing DNA sequences.
#' @param method (character, default `"dada2"`) Detection engine: `"dada2"` for
#'   de novo abundance-based detection, or `"vsearch_ref"` for reference-based
#'   UCHIME detection (requires `database`).
#' @param dada2_method (character, default `"consensus"`) Passed as `method` to
#'   `dada2::removeBimeraDenovo()` when `method = "dada2"`. One of
#'   `"consensus"`, `"pooled"`, `"per-sample"`.
#' @param database (character, default `NULL`) Path to the reference FASTA
#'   database used by vsearch. Required when `method = "vsearch_ref"`.
#' @param vsearchpath (character, default `"vsearch"`) Path to the vsearch
#'   executable.
#' @param vsearch_args (character, default `"--fasta_width 0"`) Extra arguments
#'   passed verbatim to vsearch.
#' @param keep_temporary_files (logical, default `FALSE`) Keep the temporary
#'   FASTA files written to `tempdir()` (vsearch only).
#' @param verbose (logical, default `TRUE`) Print a one-line detection summary.
#' @param ... Additional arguments passed to `dada2::removeBimeraDenovo()`.
#'
#' @return A [contam_tbl] with one row per chimeric taxon, `method = "chimera"`,
#'   and the evidence columns `chimera_method` (`"dada2"` or `"vsearch_ref"`)
#'   and `n_reads`. Empty when no chimera is detected.
#' @export
#' @author Adrien Taudière
#'
#' @seealso [filter_contam_pq()] to remove flagged taxa, [create_chimera_pq()]
#'   to generate synthetic chimeras for benchmarking.
#'
#' @examplesIf rlang::is_installed("dada2")
#' library(MiscMetabar)
#' data(data_fungi)
#' \donttest{
#' flagged <- identify_contam_chimera_pq(data_fungi, method = "dada2")
#' flagged
#' }
identify_contam_chimera_pq <- function(
  physeq,
  method = c("dada2", "vsearch_ref"),
  dada2_method = "consensus",
  database = NULL,
  vsearchpath = "vsearch",
  vsearch_args = "--fasta_width 0",
  keep_temporary_files = FALSE,
  verbose = TRUE,
  ...
) {
  method <- match.arg(method)
  MiscMetabar::verify_pq(physeq)
  if (is.null(phyloseq::refseq(physeq, errorIfNULL = FALSE))) {
    stop(
      "phyloseq object must have a refseq slot containing DNA sequences.",
      call. = FALSE
    )
  }
  taxa_nms <- phyloseq::taxa_names(physeq)

  if (method == "dada2") {
    if (!requireNamespace("dada2", quietly = TRUE)) {
      stop("Package 'dada2' is required for method = 'dada2'.", call. = FALSE)
    }
    df_seq_ab <- data.frame(
      sequence = as.character(phyloseq::refseq(physeq)),
      abundance = phyloseq::taxa_sums(physeq)
    )
    otu_nochim <- dada2::removeBimeraDenovo(
      df_seq_ab,
      method = dada2_method,
      verbose = FALSE,
      ...
    )
    kept_seqs <- rownames(otu_nochim)
    kept_taxa <- names(phyloseq::refseq(physeq))[
      names(phyloseq::refseq(physeq)) %in% kept_seqs
    ]
    chimeric_taxa <- taxa_nms[!taxa_nms %in% kept_taxa]
  } else {
    if (!requireNamespace("Biostrings", quietly = TRUE)) {
      stop("Package 'Biostrings' is required.", call. = FALSE)
    }
    if (is.null(database) || !file.exists(database)) {
      stop(
        "method = 'vsearch_ref' requires an existing `database` FASTA path.",
        call. = FALSE
      )
    }
    dna <- phyloseq::refseq(physeq)
    names(dna) <- paste0("Taxa", seq_along(dna))

    tmp_fasta <- tempfile(pattern = "chimera_ref_input_", fileext = ".fasta")
    tmp_nonchim <- tempfile(pattern = "chimera_ref_nochim_", fileext = ".fasta")
    tmp_chim <- tempfile(pattern = "chimera_ref_chim_", fileext = ".fasta")
    tmp_border <- tempfile(pattern = "chimera_ref_border_", fileext = ".fasta")
    Biostrings::writeXStringSet(dna, tmp_fasta)

    ret <- system2(
      vsearchpath,
      paste(
        "--uchime_ref",
        tmp_fasta,
        "--db",
        database,
        "--nonchimeras",
        tmp_nonchim,
        "--chimeras",
        tmp_chim,
        "--borderline",
        tmp_border,
        vsearch_args
      ),
      stdout = TRUE,
      stderr = TRUE
    )
    exit_status <- attr(ret, "status")
    if (!is.null(exit_status) && exit_status != 0) {
      stop("vsearch failed with exit code ", exit_status, call. = FALSE)
    }
    chim_seqs <- Biostrings::readDNAStringSet(tmp_chim)
    chimeric_idx <- as.integer(sub("^Taxa", "", names(chim_seqs)))
    chimeric_taxa <- taxa_nms[chimeric_idx]

    if (!keep_temporary_files) {
      unlink(c(tmp_fasta, tmp_nonchim, tmp_chim, tmp_border))
    }
  }

  chimeric_taxa <- chimeric_taxa[!is.na(chimeric_taxa)]
  chimera_engine <- method
  flagged <- tibble::tibble(
    taxon = chimeric_taxa,
    method = "chimera",
    chimera_method = chimera_engine,
    n_reads = phyloseq::taxa_sums(physeq)[chimeric_taxa]
  )

  if (verbose) {
    cli::cli_alert_info(
      "chimera ({method}): flagged {nrow(flagged)} of {length(taxa_nms)} taxa."
    )
  }

  new_contam_tbl(flagged)
}
