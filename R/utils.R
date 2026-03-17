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
