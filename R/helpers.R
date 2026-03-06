# Helper functions for phyloseq manipulation

################################################################################
#' Calculate taxa prevalence
#'
#' @description
#' Calculate the prevalence (number of samples in which a taxon is present)
#' for each taxon in a phyloseq object.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param threshold (numeric, default 0) Minimum abundance to consider a taxon
#'   as present in a sample. Values > threshold are counted as present.
#'
#' @return A named numeric vector (class integer) with prevalence for each taxon.
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' library(MiscMetabar)
#' # Get prevalence for each taxon
#' prev <- taxa_prevalence(data_fungi)
#' head(prev)
#'
#' # Get prevalence with minimum abundance of 10
#' prev10 <- taxa_prevalence(data_fungi, threshold = 10)
#' head(prev10)
#'
#' # Use in filter_taxa_pq
#' # Keep taxa present in at least 5 samples
#' filter_taxa_pq(data_fungi, taxa_prevalence(.) >= 5)
taxa_prevalence <- function(physeq, threshold = 0) {
  otu <- MiscMetabar::taxa_as_rows(physeq)@otu_table
  prev <- as.integer(rowSums(otu > threshold))
  names(prev) <- phyloseq::taxa_names(physeq)
  return(prev)
}

################################################################################
#' Decontaminate based on negative/blank control samples
#'
#' @description
#' Remove potential contaminants by setting OTU values to 0 when they are at or
#' below the level observed in negative/blank control samples for that particular OTU.
#' If multiple control are available, for each taxon a threshold is computed from the
#' control samples using a summary function (default: `max`). Occurrences in
#' non-control samples that are at or below this threshold are set to 0.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param control_condition <data-masking> An expression evaluated on sample_data
#'   that returns TRUE for control samples. Use `.` to refer to the phyloseq
#'   object (e.g., `sample_type == "negative"`, `is_control == TRUE`).
#' @param fun (function, default `max`) A function to summarize the control
#'   sample values for each taxon (or globally if `global_threshold = TRUE`).
#'   Common choices: `max` (most conservative, default), `mean`, `median`, or a
#'   custom function. The function should take a numeric vector and return a
#'   single value.
#' @param global_threshold (logical, default FALSE) If TRUE, compute a single
#'   global threshold from all control occurrences instead of per-taxon thresholds.
#'   This applies `fun` to all values in the control samples.
#' @param remove_controls (logical, default FALSE) Whether to remove the control
#'   samples from the output phyloseq object after decontamination.
#' @param clean_phyloseq_object (logical, default TRUE) Whether to clean the
#'   resulting phyloseq object using `clean_pq()` to remove empty taxa/samples.
#' @param verbose (logical, default TRUE) Whether to print additional information.
#'
#' @return A phyloseq object with decontaminated OTU values.
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' library(MiscMetabar)
#' # Add a mock control column for demonstration using the 3 samples with lowest
#' # total abundance as controls
#' pq <- mutate_samdata_pq(data_fungi, is_control = sample_sums(.) < sort(sample_sums(.))[3])
#'
#' # Decontaminate using max of controls as threshold (per-taxon)
#' decontam_sam_control(pq, is_control)
#'
#' # Use a global threshold (single value for all taxa)
#' decontam_sam_control(pq, is_control, global_threshold = TRUE)
#'
#' # Use mean instead of max (less conservative)
#' decontam_sam_control(pq, is_control, fun = mean)
#'
#' # Keep control samples in output
#' decontam_sam_control(pq, is_control, remove_controls = FALSE)
#'
#' # Use a custom function (e.g., 2x the max)
#' decontam_sam_control(pq, is_control, fun = \(x) 2 * max(x))
decontam_sam_control <- function(physeq,
                                 control_condition,
                                 fun = max,
                                 global_threshold = FALSE,
                                 remove_controls = FALSE,
                                 clean_phyloseq_object = TRUE,
                                 verbose = TRUE) {
  MiscMetabar::verify_pq(physeq)

  # Identify control samples using data masking
  mask <- build_sample_data_mask(physeq)
  condition_quo <- rlang::enquo(control_condition)
  is_control <- rlang::eval_tidy(condition_quo, data = mask)

  if (!is.logical(is_control)) {
    stop("control_condition must evaluate to a logical vector.")
  }
  if (sum(is_control) == 0) {
    stop("No samples match the control_condition.")
  }

  control_samples <- phyloseq::sample_names(physeq)[is_control]
  non_control_samples <- phyloseq::sample_names(physeq)[!is_control]

  if (length(non_control_samples) == 0) {
    stop("All samples are controls. No samples to decontaminate.")
  }

  # Get OTU matrix with taxa as rows
  tar <- phyloseq::taxa_are_rows(physeq)
  otu <- as(phyloseq::otu_table(physeq), "matrix")
  if (!tar) otu <- t(otu)

  # Calculate threshold from control samples
  control_otu <- otu[, control_samples, drop = FALSE]

  if (global_threshold) {
    # Single global threshold from all control values
    threshold_value <- fun(as.vector(control_otu))
    threshold_matrix <- matrix(threshold_value, nrow = nrow(otu), ncol = ncol(otu))
  } else {
    # Per-taxon thresholds
    thresholds <- apply(control_otu, 1, fun)
    threshold_matrix <- matrix(thresholds, nrow = nrow(otu), ncol = ncol(otu), byrow = FALSE)
  }

  # Set values to 0 where value <= threshold (only for non-control samples)
  non_control_cols <- colnames(otu) %in% non_control_samples
  otu[, non_control_cols][otu[, non_control_cols] <= threshold_matrix[, non_control_cols]] <- 0

  # Restore original orientation
  if (!tar) otu <- t(otu)

  # Update phyloseq object
  new_physeq <- physeq
  new_physeq@otu_table <- phyloseq::otu_table(otu, taxa_are_rows = tar)

  # Remove control samples if requested
  if (remove_controls) {
    new_physeq <- phyloseq::prune_samples(non_control_samples, new_physeq)
  }
  if (clean_phyloseq_object) {
    new_physeq <- MiscMetabar::clean_pq(new_physeq, silent = TRUE)
  }

  if (verbose) {
    threshold_type <- if (global_threshold) "global" else "per-taxon"
    threshold_info <- if (global_threshold) {
      paste0("\nGlobal threshold value: ", fun(as.vector(control_otu)))
    } else {
      ""
    }
    message(
      "Decontamination complete.",
      "\nThreshold type: ", threshold_type,
      "\nNumber of control samples: ", length(control_samples),
      "\nNumber of non-control samples: ", length(non_control_samples),
      "\nFunction used for threshold: ", deparse(substitute(fun)),
      threshold_info,
      "\nRemove control samples: ", remove_controls,
      "\nNumber of sequences before: ", sum(phyloseq::otu_table(physeq)),
      "\nNumber of sequences after: ", sum(phyloseq::otu_table(new_physeq)),
      "\nNumber of sequences discarded: ",
        sum(phyloseq::otu_table(physeq)) - sum(phyloseq::otu_table(new_physeq)),
      "\nNumber of occurrences discarded: ",
        sum(phyloseq::otu_table(physeq) > 0) - sum(phyloseq::otu_table(new_physeq) > 0),
      "\nNumber of taxa before: ", phyloseq::ntaxa(physeq),
      "\nNumber of taxa after: ", phyloseq::ntaxa(new_physeq),
      "\nNumber of discarded taxa: ", phyloseq::ntaxa(physeq) - phyloseq::ntaxa(new_physeq)
    )
  }

  return(new_physeq)
}

################################################################################
#' Decontaminate based on control taxa
#'
#' @description
#' Remove potential contaminants by using known control taxa (e.g., spike-ins,
#' synthetic sequences) to estimate background contamination levels. For each
#' sample, a threshold is computed from the control taxa using a summary function
#' (default: `max`). Occurrences of non-control taxa that are at or below this
#' threshold are set to 0.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param control_condition <data-masking> An expression evaluated on tax_table
#'   that returns TRUE for control taxa. Use `.` to refer to the phyloseq object.
#'   Examples: `Genus == "Tintelnotia"`, `Family == "Mitochondria"`,
#'   `taxa_names(.) %in% c("ASV1", "ASV2")`.
#' @param fun (function, default `max`) A function to summarize the control
#'   taxa values for each sample (or globally if `global_threshold = TRUE`).
#'   Common choices: `max` (most conservative, default), `mean`, `median`.
#' @param global_threshold (logical, default FALSE) If TRUE, compute a single
#'   global threshold from all control taxa occurrences instead of per-sample
#'   thresholds.
#' @param remove_control_taxa (logical, default TRUE) Whether to remove the
#'   control taxa from the output phyloseq object after decontamination.
#' @param clean_phyloseq_object (logical, default TRUE) Whether to clean the
#'   resulting phyloseq object using `clean_pq()` to remove empty taxa/samples.
#' @param verbose (logical, default TRUE) Whether to print additional information.
#'
#' @return A phyloseq object with decontaminated OTU values.
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' library(MiscMetabar)
#' # Using a condition on tax_table (e.g., select by Genus)
#' decontam_taxa_control(data_fungi, Genus == "Tintelnotia")
#'
#' # Using taxa names directly
#' control_taxa <- phyloseq::taxa_names(data_fungi)[1:2]
#' decontam_taxa_control(data_fungi, taxa_names(.) %in% control_taxa)
#'
#' # Use a global threshold
#' decontam_taxa_control(data_fungi, Genus == "Tintelnotia", global_threshold = TRUE)
#'
#' # Keep control taxa in output
#' decontam_taxa_control(data_fungi, Genus == "Tintelnotia", remove_control_taxa = FALSE)
decontam_taxa_control <- function(physeq,
                                  control_condition,
                                  fun = max,
                                  global_threshold = FALSE,
                                  remove_control_taxa = TRUE,
                                  clean_phyloseq_object = TRUE,
                                  verbose = TRUE) {
  MiscMetabar::verify_pq(physeq)

  all_taxa <- phyloseq::taxa_names(physeq)

  # Identify control taxa using data masking on tax_table
  mask <- build_taxa_data_mask(physeq)
  condition_quo <- rlang::enquo(control_condition)
  is_control <- rlang::eval_tidy(condition_quo, data = mask)

  if (!is.logical(is_control)) {
    stop("control_condition must evaluate to a logical vector.")
  }
  if (sum(is_control, na.rm = TRUE) == 0) {
    stop("No taxa match the control_condition.")
  }

  # Handle NA as FALSE

  is_control[is.na(is_control)] <- FALSE

  control_taxa_names <- all_taxa[is_control]
  non_control_taxa <- all_taxa[!is_control]

  if (length(non_control_taxa) == 0) {
    stop("All taxa are control taxa. No taxa to decontaminate.")
  }

  # Get OTU matrix with taxa as rows
  tar <- phyloseq::taxa_are_rows(physeq)
  otu <- as(phyloseq::otu_table(physeq), "matrix")
  if (!tar) otu <- t(otu)

  # Calculate threshold from control taxa
  control_otu <- otu[control_taxa_names, , drop = FALSE]

  if (global_threshold) {
    # Single global threshold from all control taxa values
    threshold_value <- fun(as.vector(control_otu))
    threshold_matrix <- matrix(threshold_value, nrow = nrow(otu), ncol = ncol(otu))
  } else {
    # Per-sample thresholds (apply fun to each column = sample)
    thresholds <- apply(control_otu, 2, fun)
    threshold_matrix <- matrix(thresholds, nrow = nrow(otu), ncol = ncol(otu), byrow = TRUE)
  }

  # Set values to 0 where value <= threshold (only for non-control taxa)
  non_control_rows <- rownames(otu) %in% non_control_taxa
  otu[non_control_rows, ][otu[non_control_rows, ] <= threshold_matrix[non_control_rows, ]] <- 0

  # Restore original orientation
  if (!tar) otu <- t(otu)

  # Update phyloseq object
  new_physeq <- physeq
  new_physeq@otu_table <- phyloseq::otu_table(otu, taxa_are_rows = tar)

  # Remove control taxa if requested
  if (remove_control_taxa) {
    new_physeq <- phyloseq::prune_taxa(non_control_taxa, new_physeq)
  }
  if (clean_phyloseq_object) {
    new_physeq <- MiscMetabar::clean_pq(new_physeq, silent = TRUE)
  }

  if (verbose) {
    threshold_type <- if (global_threshold) "global" else "per-sample"
    threshold_info <- if (global_threshold) {
      paste0("\nGlobal threshold value: ", fun(as.vector(control_otu)))
    } else {
      ""
    }
    message(
      "Decontamination complete.",
      "\nThreshold type: ", threshold_type,
      "\nNumber of control taxa: ", length(control_taxa_names),
      "\nNumber of non-control taxa: ", length(non_control_taxa),
      "\nFunction used for threshold: ", deparse(substitute(fun)),
      threshold_info,
      "\nRemove control taxa: ", remove_control_taxa,
      "\nNumber of sequences before: ", sum(phyloseq::otu_table(physeq)),
      "\nNumber of sequences after: ", sum(phyloseq::otu_table(new_physeq)),
      "\nNumber of sequences discarded: ",
        sum(phyloseq::otu_table(physeq)) - sum(phyloseq::otu_table(new_physeq)),
      "\nNumber of occurrences discarded: ",
        sum(phyloseq::otu_table(physeq) > 0) - sum(phyloseq::otu_table(new_physeq) > 0),
      "\nNumber of taxa before: ", phyloseq::ntaxa(physeq),
      "\nNumber of taxa after: ", phyloseq::ntaxa(new_physeq),
      "\nNumber of discarded taxa: ", phyloseq::ntaxa(physeq) - phyloseq::ntaxa(new_physeq)
    )
  }

  return(new_physeq)
}
