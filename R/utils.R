# Internal utility functions for tidypq
# These functions are not exported

#' Build a data mask for sample-level operations
#'
#' Creates a data mask containing sample_data columns and a `.` pronoun
#' that refers to the phyloseq object for use with phyloseq functions.
#'
#' @param physeq A phyloseq object.
#' @return A data mask for use with rlang::eval_tidy.
#' @keywords internal
#' @noRd
build_sample_data_mask <- function(physeq) {
  sam_df <- as.data.frame(phyloseq::sample_data(physeq))
  data_env <- rlang::new_environment(as.list(sam_df))
  data_env[["."]] <- physeq
  rlang::new_data_mask(data_env)
}

#' Build a data mask for taxa-level operations
#'
#' Creates a data mask containing tax_table columns and a `.` pronoun
#' that refers to the phyloseq object for use with phyloseq functions.
#'
#' @param physeq A phyloseq object.
#' @return A data mask for use with rlang::eval_tidy.
#' @keywords internal
#' @noRd
build_taxa_data_mask <- function(physeq) {
  tax_df <- as.data.frame(phyloseq::tax_table(physeq))
  data_env <- rlang::new_environment(as.list(tax_df))
  data_env[["."]] <- physeq
  rlang::new_data_mask(data_env)
}

#' Evaluate multiple quosures and combine with logical AND
#'
#' @param quos A list of quosures.
#' @param mask A data mask.
#' @return A logical vector.
#' @keywords internal
#' @noRd
eval_filter_quos <- function(quos, mask) {
  if (length(quos) == 0) {
    return(TRUE)
  }
  results <- lapply(quos, rlang::eval_tidy, data = mask)
  Reduce(`&`, results)
}

#' Build summary matrices from an OTU matrix (taxa as rows)
#'
#' Returns a named list of matrices matching OTU dimensions with per-sample
#' and per-taxon summary statistics. Used by occurrence-level functions.
#'
#' @param otu A numeric matrix with taxa as rows and samples as columns.
#' @param include_median If TRUE, also compute sample_median and taxon_median.
#' @return A named list of matrices.
#' @keywords internal
#' @noRd
build_otu_summary_matrices <- function(otu, include_median = FALSE) {
  nr <- nrow(otu)
  nc <- ncol(otu)

  out <- list(
    `.` = otu,
    sample_total = matrix(colSums(otu), nrow = nr, ncol = nc, byrow = TRUE),
    taxon_total = matrix(rowSums(otu), nrow = nr, ncol = nc, byrow = FALSE),
    sample_mean = matrix(colMeans(otu), nrow = nr, ncol = nc, byrow = TRUE),
    taxon_mean = matrix(rowMeans(otu), nrow = nr, ncol = nc, byrow = FALSE)
  )

  if (include_median) {
    out$sample_median <- matrix(
      apply(otu, 2, stats::median),
      nrow = nr,
      ncol = nc,
      byrow = TRUE
    )
    out$taxon_median <- matrix(
      apply(otu, 1, stats::median),
      nrow = nr,
      ncol = nc,
      byrow = FALSE
    )
  }

  out
}

#' Print a standardized decontamination summary message
#'
#' @param physeq Original phyloseq object.
#' @param new_physeq Decontaminated phyloseq object.
#' @param threshold_type Character, e.g. "global" or "per-taxon".
#' @param n_control Number of control units.
#' @param n_non_control Number of non-control units.
#' @param fun_name Character name of the summary function.
#' @param global_threshold_value Numeric or NULL.
#' @param control_label "samples" or "taxa".
#' @param removed Logical, whether controls were removed.
#' @keywords internal
#' @noRd
decontam_message <- function(
  physeq,
  new_physeq,
  threshold_type,
  n_control,
  n_non_control,
  fun_name,
  global_threshold_value = NULL,
  control_label,
  removed
) {
  seqs_before <- sum(phyloseq::otu_table(physeq))
  seqs_after <- sum(phyloseq::otu_table(new_physeq))
  occ_before <- sum(phyloseq::otu_table(physeq) > 0)
  occ_after <- sum(phyloseq::otu_table(new_physeq) > 0)
  taxa_before <- phyloseq::ntaxa(physeq)
  taxa_after <- phyloseq::ntaxa(new_physeq)

  threshold_info <- if (!is.null(global_threshold_value)) {
    paste0(
      "\nGlobal threshold value: ",
      global_threshold_value
    )
  } else {
    ""
  }

  message(
    "Decontamination complete.",
    "\nThreshold type: ",
    threshold_type,
    "\nControl ",
    control_label,
    ": ",
    n_control,
    "\nNon-control ",
    control_label,
    ": ",
    n_non_control,
    "\nThreshold function: ",
    fun_name,
    threshold_info,
    "\nControls removed: ",
    removed,
    "\nSequences: ",
    seqs_before,
    " -> ",
    seqs_after,
    " (-",
    seqs_before - seqs_after,
    ")",
    "\nOccurrences: ",
    occ_before,
    " -> ",
    occ_after,
    " (-",
    occ_before - occ_after,
    ")",
    "\nTaxa: ",
    taxa_before,
    " -> ",
    taxa_after,
    " (-",
    taxa_before - taxa_after,
    ")"
  )
}
