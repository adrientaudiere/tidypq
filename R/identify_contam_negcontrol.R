# Detector: classify negative-control taxa into contamination sub-types

################################################################################
#' Detect contaminants from negative-control occurrence patterns
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Examines taxa detected in negative-control samples and assigns each a
#' contamination sub-type from read abundance and occurrence patterns:
#'
#' * **artifact**: very low total reads AND present in very few samples --
#'   likely sequencing noise or index hopping.
#' * **lab_contaminant**: predominantly found in negative controls relative to
#'   real samples -- likely introduced during library preparation.
#' * **sample_contaminant**: detected in controls but also widespread in real
#'   samples -- may represent genuine biological signal.
#'
#' This is a pure detector: it flags taxa and returns a [contam_tbl] without
#' modifying `physeq`. Because `sample_contaminant` may be genuine signal, only
#' `artifact` and `lab_contaminant` are flagged by default (see
#' `flag_categories`). Use [filter_contam_pq()] to remove the flagged taxa.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param neg_control <data-masking> An expression evaluated on `sample_data`
#'   returning `TRUE` for negative-control samples (e.g. `is_control == TRUE`).
#'   Use `.` to refer to the phyloseq object.
#' @param min_reads_artifact (numeric, default `10`) Maximum total reads for a
#'   taxon to be considered an artifact.
#' @param max_samples_artifact (integer, default `2`) Maximum total number of
#'   samples in which an artifact taxon may appear.
#' @param max_ratio_lab_contam (numeric, default `0.2`) Maximum ratio of
#'   non-control to control sample occurrences for the lab-contaminant category.
#' @param min_neg_samples_lab (integer, default `2`) Minimum number of control
#'   samples a lab-contaminant taxon must appear in.
#' @param flag_categories (character, default `c("artifact",
#'   "lab_contaminant")`) Which sub-types to include in the returned
#'   `contam_tbl`. Add `"sample_contaminant"` to flag those as well.
#' @param verbose (logical, default `TRUE`) Print a one-line detection summary.
#'
#' @return A [contam_tbl] with one row per flagged taxon, `method =
#'   "negcontrol"`, and the evidence columns `subtype` (`artifact`,
#'   `lab_contaminant`, or `sample_contaminant`), `total_reads`, `reads_in_neg`,
#'   `reads_in_samples`, `n_neg_samples`, `n_non_neg_samples`, and
#'   `ratio_non_neg_to_neg`. Empty when no eligible taxon is found.
#' @export
#' @author Adrien Taudière
#'
#' @seealso [filter_contam_pq()] to remove flagged taxa.
#'
#' @examples
#' library(MiscMetabar)
#' data(data_fungi)
#' \donttest{
#' pq <- mutate_samdata_pq(
#'   data_fungi,
#'   is_control = sample_sums(.) < sort(sample_sums(.))[4]
#' )
#' flagged <- identify_contam_negcontrol_pq(pq, is_control)
#' flagged
#' }
identify_contam_negcontrol_pq <- function(
  physeq,
  neg_control,
  min_reads_artifact = 10,
  max_samples_artifact = 2,
  max_ratio_lab_contam = 0.2,
  min_neg_samples_lab = 2,
  flag_categories = c("artifact", "lab_contaminant"),
  verbose = TRUE
) {
  MiscMetabar::verify_pq(physeq)
  flag_categories <- match.arg(
    flag_categories,
    choices = c("artifact", "lab_contaminant", "sample_contaminant"),
    several.ok = TRUE
  )

  mask <- build_sample_data_mask(physeq)
  nc_quo <- rlang::enquo(neg_control)
  nc_vec <- rlang::eval_tidy(nc_quo, data = mask)
  if (!is.logical(nc_vec)) {
    cli::cli_abort("{.arg neg_control} must evaluate to a logical vector.")
  }
  nc_vec[is.na(nc_vec)] <- FALSE
  if (sum(nc_vec) == 0) {
    cli::cli_abort("No samples match the {.arg neg_control} expression.")
  }

  otu <- as(phyloseq::otu_table(physeq), "matrix")
  if (!phyloseq::taxa_are_rows(physeq)) {
    otu <- t(otu)
  }
  storage.mode(otu) <- "double"

  otu_neg <- otu[, which(nc_vec), drop = FALSE]
  otu_non_neg <- otu[, which(!nc_vec), drop = FALSE]
  taxa_in_neg <- rownames(otu_neg)[rowSums(otu_neg > 0) >= 1]

  if (length(taxa_in_neg) == 0) {
    if (verbose) {
      cli::cli_alert_info("negcontrol: no taxa found in negative controls.")
    }
    return(new_contam_tbl(tibble::tibble(
      taxon = character(),
      method = character(),
      subtype = character()
    )))
  }

  clf <- tibble::tibble(
    taxon = taxa_in_neg,
    method = "negcontrol",
    subtype = NA_character_,
    total_reads = rowSums(otu[taxa_in_neg, , drop = FALSE]),
    reads_in_neg = rowSums(otu_neg[taxa_in_neg, , drop = FALSE]),
    reads_in_samples = rowSums(otu_non_neg[taxa_in_neg, , drop = FALSE]),
    n_neg_samples = rowSums(otu_neg[taxa_in_neg, , drop = FALSE] > 0),
    n_non_neg_samples = rowSums(otu_non_neg[taxa_in_neg, , drop = FALSE] > 0)
  )
  clf$ratio_non_neg_to_neg <- clf$n_non_neg_samples / clf$n_neg_samples

  is_artifact <- clf$total_reads <= min_reads_artifact &
    (clf$n_neg_samples + clf$n_non_neg_samples) <= max_samples_artifact
  clf$subtype[is_artifact] <- "artifact"

  is_lab <- is.na(clf$subtype) &
    clf$ratio_non_neg_to_neg <= max_ratio_lab_contam &
    clf$n_neg_samples >= min_neg_samples_lab
  clf$subtype[is_lab] <- "lab_contaminant"

  clf$subtype[is.na(clf$subtype)] <- "sample_contaminant"

  flagged <- clf[clf$subtype %in% flag_categories, , drop = FALSE]
  flagged <- flagged[order(flagged$subtype, -flagged$total_reads), ]

  if (verbose) {
    counts <- table(factor(flagged$subtype, levels = flag_categories))
    cli::cli_alert_info(
      "negcontrol: flagged {nrow(flagged)} taxa ({paste(names(counts), counts, sep = '=', collapse = ', ')})."
    )
  }

  new_contam_tbl(flagged)
}
