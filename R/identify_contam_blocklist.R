# Detector: taxa whose genus belongs to a known-contaminant blocklist

################################################################################
#' Detect taxa whose genus belongs to a known-contaminant blocklist
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Flags taxa whose genus matches a blocklist of known reagent and laboratory
#' contaminants. A taxonomy-only method requiring no negative controls, best
#' used as a complementary, curation-oriented filter -- review flagged taxa
#' before removal, as several blocklist genera also occur as genuine
#' environmental taxa.
#'
#' This is a pure detector: it flags taxa and returns a [contam_tbl] without
#' modifying `physeq`. Use [filter_contam_pq()] to remove the flagged taxa.
#'
#' @param physeq (phyloseq, required) A phyloseq object with a `tax_table`.
#' @param genus_blocklist (character, default [known_contaminant_genera()])
#'   Vector of genus names to flag.
#' @param tax_rank (character, default `"Genus"`) Name of the `tax_table`
#'   column holding genus-level assignments.
#' @param extra_genera (character, default `NULL`) Additional genus names to
#'   append to `genus_blocklist`.
#' @param ignore_case (logical, default `TRUE`) If `TRUE`, matching is
#'   case-insensitive.
#' @param verbose (logical, default `TRUE`) Print a one-line detection summary.
#'
#' @return A [contam_tbl] with one row per flagged taxon, `method = "blocklist"`,
#'   and the evidence columns `genus`, `total_reads`, and `prevalence`. Empty
#'   when no genus matches the blocklist.
#' @export
#' @author Adrien Taudière
#'
#' @seealso [known_contaminant_genera()] for the default blocklist,
#'   [filter_contam_pq()] to remove flagged taxa.
#'
#' @examples
#' library(MiscMetabar)
#' data(data_fungi)
#' \donttest{
#' # data_fungi is fungal, so the default bacterial blocklist matches nothing;
#' # extend it with a genus present in the data for illustration.
#' identify_contam_blocklist_pq(data_fungi, extra_genera = "Mortierella")
#' }
identify_contam_blocklist_pq <- function(
  physeq,
  genus_blocklist = known_contaminant_genera(),
  tax_rank = "Genus",
  extra_genera = NULL,
  ignore_case = TRUE,
  verbose = TRUE
) {
  MiscMetabar::verify_pq(physeq)

  if (is.null(phyloseq::tax_table(physeq, errorIfNULL = FALSE))) {
    cli::cli_abort("{.arg physeq} must contain a {.field tax_table}.")
  }
  tax <- taxa_table_df(physeq)
  if (!tax_rank %in% colnames(tax)) {
    cli::cli_abort(c(
      "Column {.val {tax_rank}} not found in the {.field tax_table}.",
      "i" = "Available ranks: {.val {colnames(tax)}}."
    ))
  }

  blocklist <- unique(c(genus_blocklist, extra_genera))
  blocklist <- blocklist[!is.na(blocklist) & nzchar(blocklist)]

  genus <- tax[[tax_rank]]
  if (ignore_case) {
    is_contaminant <- tolower(genus) %in% tolower(blocklist)
  } else {
    is_contaminant <- genus %in% blocklist
  }
  is_contaminant[is.na(genus)] <- FALSE

  otu <- as(phyloseq::otu_table(physeq), "matrix")
  if (!phyloseq::taxa_are_rows(physeq)) {
    otu <- t(otu)
  }
  storage.mode(otu) <- "double"
  otu <- otu[phyloseq::taxa_names(physeq), , drop = FALSE]

  bad_idx <- which(is_contaminant)
  flagged <- tibble::tibble(
    taxon = phyloseq::taxa_names(physeq)[bad_idx],
    method = "blocklist",
    genus = genus[bad_idx],
    total_reads = rowSums(otu)[bad_idx],
    prevalence = rowSums(otu > 0)[bad_idx]
  )
  flagged <- flagged[order(-flagged$total_reads), , drop = FALSE]

  if (verbose) {
    n_total <- phyloseq::ntaxa(physeq)
    cli::cli_alert_info(
      "blocklist: flagged {nrow(flagged)} of {n_total} taxa from {length(blocklist)} blocklist genera."
    )
  }

  new_contam_tbl(flagged)
}
