# Convenience verb: detect and remove chimeric taxa in one call

################################################################################
#' Remove chimeric taxa in one step
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' A convenience wrapper that runs [identify_contam_chimera_pq()] and then
#' [filter_contam_pq()], returning the chimera-free phyloseq object directly.
#' Chimera removal is a near-ubiquitous cleaning step, so it gets a dedicated
#' verb; the other contaminant methods are composed explicitly via their
#' `identify_contam_*_pq()` detector and `filter_contam_pq()`.
#'
#' @param physeq (phyloseq, required) A phyloseq object with a populated
#'   `@@refseq` slot containing DNA sequences.
#' @param method (character, default `"dada2"`) Detection engine, passed to
#'   [identify_contam_chimera_pq()]: `"dada2"` (de novo) or `"vsearch_ref"`
#'   (reference-based UCHIME, requires `database`).
#' @param clean_phyloseq_object (logical, default `TRUE`) Whether to clean the
#'   result with `MiscMetabar::clean_pq()`.
#' @param verbose (logical, default `TRUE`) Whether to report detection and
#'   removal.
#' @param ... Additional arguments passed to [identify_contam_chimera_pq()]
#'   (e.g. `dada2_method`, `database`, `vsearchpath`, `vsearch_args`).
#'
#' @return A phyloseq object with chimeric taxa removed (unchanged when none are
#'   detected).
#' @export
#' @author Adrien Taudière
#'
#' @seealso [identify_contam_chimera_pq()] to inspect chimeras before removal,
#'   [filter_contam_pq()] for the general removal verb.
#'
#' @examplesIf rlang::is_installed("dada2")
#' library(MiscMetabar)
#' data(data_fungi)
#' \donttest{
#' data_nochim <- filter_chimera_pq(data_fungi, method = "dada2")
#' }
filter_chimera_pq <- function(
  physeq,
  method = c("dada2", "vsearch_ref"),
  clean_phyloseq_object = TRUE,
  verbose = TRUE,
  ...
) {
  flagged <- identify_contam_chimera_pq(
    physeq,
    method = method,
    verbose = verbose,
    ...
  )
  filter_contam_pq(
    physeq,
    flagged,
    clean_phyloseq_object = clean_phyloseq_object,
    verbose = verbose
  )
}
