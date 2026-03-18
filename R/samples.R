# Sample-scale tidyverse verbs for phyloseq objects

################################################################################
#' Filter samples in a phyloseq object
#'
#' @description
#' Filter samples using data masking on sample_data. Supports the `.` pronoun
#' to refer to the phyloseq object for use with functions like `sample_sums()`.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param ... <data-masking> Expressions that return a logical value, evaluated
#'   in the context of sample_data. Multiple conditions are combined with `&`.
#'   Use `.` to refer to the phyloseq object (e.g., `sample_sums(.) > 1000`).
#' @param clean_phyloseq_object if TRUE (default), the resulting phyloseq object
#' is cleaned using `clean_pq()` to remove empty taxa/samples.
#' @return A phyloseq object with filtered samples.
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' library(MiscMetabar)
#' # Filter by sample metadata
#' filter_samples_pq(data_fungi, Height == "Low")
#'
#' # Filter by sequencing depth
#' filter_samples_pq(data_fungi, sample_sums(.) > 1000)
#'
#' # Combine multiple conditions
#' filter_samples_pq(data_fungi, Height == "Low", sample_sums(.) > 5000)
#'
#' # Keep samples above median abundance
#' filter_samples_pq(data_fungi, sample_sums(.) > median(sample_sums(.)))
#'
#' # Keep samples above half of the average abundance
#' filter_samples_pq(data_fungi, sample_sums(.) > sum(sample_sums(.))/phyloseq::nsamples(.)/2)
filter_samples_pq <- function(physeq, ..., clean_phyloseq_object = TRUE) {
  MiscMetabar::verify_pq(physeq)

  mask <- build_sample_data_mask(physeq)
  quos <- rlang::enquos(...)
  keep <- eval_filter_quos(quos, mask)

  samples_to_keep <- phyloseq::sample_names(physeq)[keep]

  if (length(samples_to_keep) == 0) {
    warning("No samples match the filter criteria.")
  }

  new_physeq <- phyloseq::prune_samples(samples_to_keep, physeq)
  if (clean_phyloseq_object) {
    return(MiscMetabar::clean_pq(new_physeq, silent = TRUE))
  } else {
    return(new_physeq)
  }
}

################################################################################
#' Select columns from sample_data in a phyloseq object
#'
#' @description
#' Select sample_data columns using tidyselect semantics.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param ... <tidy-select> One or more unquoted expressions separated by
#'   commas. Variable names can be used as if they were positions in the data
#'   frame, so expressions like `x:y` can be used to select a range of
#'   variables.
#'
#' @return A phyloseq object with selected sample_data columns.
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' library(MiscMetabar)
#' # Select specific columns
#' select_samdata_pq(data_fungi, Height, Time)
#'
#' # Select a range of columns
#' select_samdata_pq(data_fungi, Height:Time)
#'
#' # Exclude columns
#' select_samdata_pq(data_fungi, !Sample_id)
select_samdata_pq <- function(physeq, ...) {
  MiscMetabar::verify_pq(physeq)
  new_physeq <- physeq

  sam_df <- as_tibble(phyloseq::sample_data(physeq))
  loc <- tidyselect::eval_select(rlang::expr(c(...)), sam_df)
  sam_df <- sam_df[, loc, drop = FALSE]

  new_physeq@sam_data <- phyloseq::sample_data(sam_df)
  sample_names(new_physeq) <- phyloseq::sample_names(physeq) # Preserve sample names
  MiscMetabar::verify_pq(new_physeq)
  return(new_physeq)
}

################################################################################
#' Add or modify columns in sample_data
#'
#' @description
#' Create new columns or modify existing ones in sample_data using data masking.
#' Supports the `.` pronoun to refer to the phyloseq object.
#'
#' This function only modifies the sample_data slot (columns/metadata). It cannot
#' add or remove samples. The number of samples and sample names are preserved.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param ... <data-masking> Name-value pairs. The name gives the name of the
#'   column in the output. The value must be a vector of length 1 (recycled) or
#'   exactly the same length as the number of samples. Use `.` to refer to the
#'
#'   phyloseq object.
#'
#' @return A phyloseq object with modified sample_data (same samples, modified
#'   or new columns).
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' library(MiscMetabar)
#' # Add a new column based on sequencing depth
#' mutate_samdata_pq(data_fungi, log_depth = log(sample_sums(.)))
#'
#' # Modify an existing column
#' mutate_samdata_pq(data_fungi, Height = toupper(Height))
mutate_samdata_pq <- function(physeq, ...) {
  MiscMetabar::verify_pq(physeq)
  new_physeq <- physeq

  sam_df <- as.data.frame(phyloseq::sample_data(physeq))
  original_nsamples <- nrow(sam_df)
  mask <- build_sample_data_mask(physeq)

  quos <- rlang::enquos(...)
  for (nm in names(quos)) {
    value <- rlang::eval_tidy(quos[[nm]], data = mask)
    if (length(value) != 1 && length(value) != original_nsamples) {
      stop(
        sprintf(
          "Column '%s' has length %d, but must be length 1 or %d (number of samples).",
          nm,
          length(value),
          original_nsamples
        )
      )
    }
    sam_df[[nm]] <- value
  }

  new_physeq@sam_data <- phyloseq::sample_data(sam_df)

  MiscMetabar::verify_pq(new_physeq)
  return(new_physeq)
}

################################################################################
#' Subset samples by position
#'
#' @description
#' Select samples by their integer positions.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param clean_phyloseq_object if TRUE (default), the resulting phyloseq object
#' is cleaned using `clean_pq()` to remove empty taxa/samples.
#' @param ... <data-masking> Integer row indices. Positive values select
#'   samples, negative values drop samples.
#'
#' @return A phyloseq object with selected samples.
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' library(MiscMetabar)
#' # Select first 5 samples
#' slice_samples_pq(data_fungi, 1:5)
#'
#' # Remove first 2 samples
#' slice_samples_pq(data_fungi, -(1:2))
slice_samples_pq <- function(physeq, ..., clean_phyloseq_object = TRUE) {
  MiscMetabar::verify_pq(physeq)

  sam_df <- as_tibble(phyloseq::sample_data(physeq))
  sam_df$sample_names__for_slicing <- phyloseq::sample_names(physeq)
  indices <- c(...)
  sam_df <- dplyr::slice(sam_df, indices)

  samples_to_keep <- sam_df$sample_names__for_slicing

  if (length(samples_to_keep) == 0) {
    warning("No samples selected.")
  }

  new_physeq <- phyloseq::prune_samples(samples_to_keep, physeq)
  if (clean_phyloseq_object) {
    return(MiscMetabar::clean_pq(new_physeq, silent = TRUE))
  } else {
    return(new_physeq)
  }
}

################################################################################
#' Arrange samples by column values
#'
#' @description
#' Reorder samples based on sample_data columns. Supports the `.` pronoun
#' to refer to the phyloseq object for sorting by computed values.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param clean_phyloseq_object if TRUE (default), the resulting phyloseq object
#'  is cleaned using `clean_pq()` to remove empty taxa/samples.
#' @param ... <data-masking> Variables to sort by. Use `desc()` for descending
#'   order. Use `.` to refer to the phyloseq object.
#'
#' @return A phyloseq object with reordered samples.
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' library(MiscMetabar)
#' # Arrange by a single column
#' arrange_samples_pq(data_fungi, Height)
#'
#' # Arrange by sequencing depth (descending)
#' arrange_samples_pq(data_fungi, dplyr::desc(sample_sums(.)))
arrange_samples_pq <- function(physeq, ..., clean_phyloseq_object = TRUE) {
  MiscMetabar::verify_pq(physeq)

  sam_df <- as.data.frame(phyloseq::sample_data(physeq))
  mask <- build_sample_data_mask(physeq)

  quos <- rlang::enquos(...)
  # Evaluate each expression in the mask context and add to data frame
  for (i in seq_along(quos)) {
    col_name <- paste0(".sort_", i)
    sam_df[[col_name]] <- rlang::eval_tidy(quos[[i]], data = mask)
  }

  # Sort by the temporary columns
  sort_cols <- paste0(".sort_", seq_along(quos))
  sam_df <- sam_df[do.call(order, sam_df[sort_cols]), , drop = FALSE]

  # Remove temporary columns

  sam_df <- sam_df[, !names(sam_df) %in% sort_cols, drop = FALSE]

  samples_order <- rownames(sam_df)

  new_physeq <- phyloseq::prune_samples(samples_order, physeq)
  if (clean_phyloseq_object) {
    return(MiscMetabar::clean_pq(new_physeq, silent = TRUE))
  } else {
    return(new_physeq)
  }
}

################################################################################
#' Rename columns in sample_data
#'
#' @description
#' Rename columns in sample_data using tidyselect semantics.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param ... <tidy-select> Name-value pairs where the name is the new name
#'   and the value is the old name. Use `new_name = old_name` syntax.
#'
#' @return A phyloseq object with renamed sample_data columns.
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' library(MiscMetabar)
#' # Rename a single column
#' rename_samples_pq(data_fungi, sample_height = Height)
#'
#' # Rename multiple columns
#' rename_samples_pq(data_fungi, sample_height = Height, sample_time = Time)
rename_samples_pq <- function(physeq, ...) {
  MiscMetabar::verify_pq(physeq)
  new_physeq <- physeq

  sam_df <- as_tibble(phyloseq::sample_data(physeq))
  sam_df <- dplyr::rename(sam_df, ...)

  new_physeq@sam_data <- phyloseq::sample_data(sam_df)
  phyloseq::sample_names(new_physeq) <- phyloseq::sample_names(physeq)

  MiscMetabar::verify_pq(new_physeq)
  return(new_physeq)
}
