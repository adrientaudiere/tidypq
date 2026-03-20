# Occurrence-scale tidyverse verbs for phyloseq objects

################################################################################
#' Filter occurrences in the OTU table
#'
#' @description
#' Set OTU table values to 0 based on a condition. This is useful for removing
#' singletons, low-abundance values, or other filtering operations at the
#' cell level.
#'
#' The condition is evaluated vectorized across the entire OTU matrix with
#' access to special variables (all as matrices matching OTU dimensions):
#' - `.` = cell values (the OTU matrix)
#' - `sample_total` = sum of each sample (repeated per column)
#' - `taxon_total` = sum of each taxon (repeated per row)
#' - `sample_mean` = mean of each sample (repeated per column)
#' - `taxon_mean` = mean of each taxon (repeated per row)
#'
#' Values that do not satisfy the condition are set to 0.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param condition An expression evaluated on the OTU matrix. Values where the
#'   condition is FALSE (or NA) are set to 0. Use `.` to refer to cell values.
#' @param clean_phyloseq_object if TRUE (default), the resulting phyloseq object
#'  is cleaned using `clean_pq()` to remove empty taxa/samples.
#'
#' @return A phyloseq object with filtered OTU values.
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' library(MiscMetabar)
#' # Remove singletons (abundance = 1)
#' filter_occurrences_pq(data_fungi, . > 1)
#'
#' # Keep only values above 0.01% of sample total
#' filter_occurrences_pq(data_fungi, . / sample_total > 0.0001)
#'
#' # Keep only values above taxon mean
#' filter_occurrences_pq(data_fungi, . > taxon_mean)
filter_occurrences_pq <- function(
  physeq,
  condition,
  clean_phyloseq_object = TRUE
) {
  MiscMetabar::verify_pq(physeq)

  tar <- phyloseq::taxa_are_rows(physeq)
  otu <- as(phyloseq::otu_table(physeq), "matrix")
  if (!tar) {
    otu <- t(otu)
  }

  summaries <- build_otu_summary_matrices(otu)

  # Evaluate condition vectorized
  condition_quo <- rlang::enquo(condition)
  mask_env <- rlang::new_environment(summaries)
  mask <- rlang::new_data_mask(mask_env)
  keep_matrix <- rlang::eval_tidy(condition_quo, data = mask)

  # Handle NA as FALSE
  keep_matrix[is.na(keep_matrix)] <- FALSE

  # Set values to 0 where condition is not met
  otu[!keep_matrix] <- 0

  # Restore original orientation
  if (!tar) {
    otu <- t(otu)
  }

  new_physeq <- physeq
  new_physeq@otu_table <- phyloseq::otu_table(otu, taxa_are_rows = tar)

  if (clean_phyloseq_object) {
    return(MiscMetabar::clean_pq(new_physeq, silent = TRUE))
  } else {
    return(new_physeq)
  }
}

################################################################################
#' Transform OTU table values
#'
#' @description
#' Apply a transformation to all values in the OTU table. This is useful for
#' computing relative abundances, log transformations, or other value-level
#' operations.
#'
#' The expression is evaluated vectorized across the entire OTU matrix with
#' access to special variables (all as matrices matching OTU dimensions):
#' - `.` = cell values (the OTU matrix)
#' - `sample_total` = sum of each sample (repeated per column)
#' - `taxon_total` = sum of each taxon (repeated per row)
#' - `sample_mean` = mean of each sample (repeated per column)
#' - `taxon_mean` = mean of each taxon (repeated per row)
#' - `sample_median` = median of each sample (repeated per column)
#' - `taxon_median` = median of each taxon (repeated per row)
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param expr An expression evaluated on the OTU matrix. The result replaces
#'   the original values. Use `.` to refer to cell values.
#'
#' @return A phyloseq object with transformed OTU values.
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' library(MiscMetabar)
#' # Convert to relative abundance (proportion of sample total)
#' mutate_occurrences_pq(data_fungi, . / sample_total)
#'
#' # Log transformation (adding pseudocount)
#' mutate_occurrences_pq(data_fungi, log1p(.))
#'
#' # Center by taxon mean
#' mutate_occurrences_pq(data_fungi, . - taxon_mean)
mutate_occurrences_pq <- function(physeq, expr) {
  MiscMetabar::verify_pq(physeq)

  tar <- phyloseq::taxa_are_rows(physeq)
  otu <- as(phyloseq::otu_table(physeq), "matrix")
  if (!tar) {
    otu <- t(otu)
  }

  summaries <- build_otu_summary_matrices(otu, include_median = TRUE)

  # Evaluate expression vectorized
  expr_quo <- rlang::enquo(expr)
  mask_env <- rlang::new_environment(summaries)
  mask <- rlang::new_data_mask(mask_env)
  new_otu <- rlang::eval_tidy(expr_quo, data = mask)

  # Preserve dimnames
  rownames(new_otu) <- rownames(otu)
  colnames(new_otu) <- colnames(otu)

  # Restore original orientation
  if (!tar) {
    new_otu <- t(new_otu)
  }

  new_physeq <- physeq
  new_physeq@otu_table <- phyloseq::otu_table(new_otu, taxa_are_rows = tar)

  MiscMetabar::verify_pq(new_physeq)
  new_physeq
}
