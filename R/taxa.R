# Taxa-scale tidyverse verbs for phyloseq objects

################################################################################
#' Filter taxa in a phyloseq object
#'
#' @description
#' Filter taxa using data masking on tax_table. Supports the `.` pronoun
#' to refer to the phyloseq object for use with functions like `taxa_sums()`.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param clean_phyloseq_object if TRUE (default), the resulting phyloseq object
#' is cleaned using `clean_pq()` to remove empty taxa/samples.
#' @param ... <data-masking> Expressions that return a logical value, evaluated
#'   in the context of tax_table. Multiple conditions are combined with `&`.
#'   Use `.` to refer to the phyloseq object (e.g., `taxa_sums(.) > 100`).
#'
#' @return A phyloseq object with filtered taxa.
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' library(MiscMetabar)
#' # Filter by taxonomy
#' filter_taxa_pq(data_fungi, Phylum == "Basidiomycota")
#'
#' # Filter by total abundance
#' filter_taxa_pq(data_fungi, taxa_sums(.) > 100)
#'
#' # Combine multiple conditions
#' filter_taxa_pq(data_fungi, Phylum == "Basidiomycota", taxa_sums(.) > 100)
#' 
#' # Keep taxa above median abundance
#' filter_taxa_pq(data_fungi, taxa_sums(.) > median(taxa_sums(.)))
#' 
filter_taxa_pq <- function(physeq, ..., clean_phyloseq_object = TRUE) {
  MiscMetabar::verify_pq(physeq)

  mask <- build_taxa_data_mask(physeq)
  quos <- rlang::enquos(...)
  keep <- eval_filter_quos(quos, mask)

  taxa_to_keep <- phyloseq::taxa_names(physeq)[keep]

  if (length(taxa_to_keep) == 0) {
    warning("No taxa match the filter criteria.")
  }


  new_physeq <- phyloseq::prune_taxa(taxa_to_keep, physeq)
  if (clean_phyloseq_object) {
    return(MiscMetabar::clean_pq(new_physeq, silent = TRUE))
  } else {
    return(new_physeq)
  }
}

################################################################################
#' Select columns from tax_table in a phyloseq object
#'
#' @description
#' Select tax_table columns (taxonomic ranks) using tidyselect semantics.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param ... <tidy-select> One or more unquoted expressions separated by
#'   commas. Variable names can be used as if they were positions in the data
#'   frame, so expressions like `Kingdom:Genus` can be used to select a range
#'   of taxonomic ranks.
#'
#' @return A phyloseq object with selected tax_table columns.
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' library(MiscMetabar)
#' # Select specific ranks
#' select_taxa_pq(data_fungi, Kingdom, Phylum, Class)
#'
#' # Select a range of ranks
#' select_taxa_pq(data_fungi, Kingdom:Genus)
#'
#' # Exclude ranks
#' select_taxa_pq(data_fungi, !Species)
select_taxa_pq <- function(physeq, ...) {
  MiscMetabar::verify_pq(physeq)
  new_physeq <- physeq

  tax_df <- as.data.frame(phyloseq::tax_table(physeq))
  loc <- tidyselect::eval_select(rlang::expr(c(...)), tax_df)
  tax_df <- tax_df[, loc, drop = FALSE]

  new_physeq@tax_table <- phyloseq::tax_table(as.matrix(tax_df))
  phyloseq::taxa_names(new_physeq@tax_table) <- phyloseq::taxa_names(physeq)

  MiscMetabar::verify_pq(new_physeq)
  new_physeq
}

################################################################################
#' Add or modify columns in tax_table
#'
#' @description
#' Create new columns or modify existing ones in tax_table using data masking.
#' Supports the `.` pronoun to refer to the phyloseq object.
#'
#' This function only modifies the tax_table slot (columns/taxonomic ranks). It
#' cannot add or remove taxa. The number of taxa and taxa names are preserved.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param ... <data-masking> Name-value pairs. The name gives the name of the
#'   column in the output. The value must be a vector of length 1 (recycled) or
#'   exactly the same length as the number of taxa. Use `.` to refer to the
#'   phyloseq object.
#'
#' @return A phyloseq object with modified tax_table (same taxa, modified or
#'   new columns).
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' library(MiscMetabar)
#' # Replace NA values in a column
#' mutate_taxa_pq(data_fungi, Genus = ifelse(is.na(Genus), "Unknown", Genus))
#'
#' # Add a new column based on abundance
#' mutate_taxa_pq(data_fungi, total_abundance = taxa_sums(.))
mutate_taxa_pq <- function(physeq, ...) {
  MiscMetabar::verify_pq(physeq)
  new_physeq <- physeq

  tax_df <- as.data.frame(phyloseq::tax_table(physeq))
  original_taxa <- rownames(tax_df)
  original_ntaxa <- nrow(tax_df)
  mask <- build_taxa_data_mask(physeq)

  quos <- rlang::enquos(...)
  for (nm in names(quos)) {
    value <- rlang::eval_tidy(quos[[nm]], data = mask)
    if (length(value) != 1 && length(value) != original_ntaxa) {
      stop(
        sprintf(
          "Column '%s' has length %d, but must be length 1 or %d (number of taxa).",
          nm, length(value), original_ntaxa
        )
      )
    }
    tax_df[[nm]] <- value
  }

  # Verify taxa are preserved
  if (nrow(tax_df) != original_ntaxa || !identical(rownames(tax_df), original_taxa)) {
    stop("mutate_taxa_pq cannot add or remove taxa. Use filter_taxa_pq or slice_taxa_pq instead.")
  }

  new_physeq@tax_table <- phyloseq::tax_table(as.matrix(tax_df))
  phyloseq::taxa_names(new_physeq@tax_table) <- phyloseq::taxa_names(physeq)

  MiscMetabar::verify_pq(new_physeq)
  new_physeq
}

################################################################################
#' Subset taxa by position
#'
#' @description
#' Select taxa by their integer positions.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param ... <data-masking> Integer row indices. Positive values select
#'   taxa, negative values drop taxa.
#' @param clean_phyloseq_object if TRUE (default), the resulting phyloseq object
#' is cleaned using `clean_pq()` to remove empty taxa/samples.
#' @return A phyloseq object with selected taxa.
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' library(MiscMetabar)
#' # Select first 10 taxa
#' slice_taxa_pq(data_fungi, 1:10)
#'
#' # Remove first 5 taxa
#' slice_taxa_pq(data_fungi, -(1:5))
slice_taxa_pq <- function(physeq, ..., clean_phyloseq_object = TRUE) {
  MiscMetabar::verify_pq(physeq)

  tax_df <- as.data.frame(phyloseq::tax_table(physeq))
  indices <- c(...)
  tax_df <- dplyr::slice(tax_df, indices)

  taxa_to_keep <- rownames(tax_df)

  if (length(taxa_to_keep) == 0) {
    warning("No taxa selected.")
  }

  new_physeq <- phyloseq::prune_taxa(taxa_to_keep, physeq)
  if (clean_phyloseq_object) {
    return(MiscMetabar::clean_pq(new_physeq, silent = TRUE))
  } else {
    return(new_physeq)
  }
}

################################################################################
#' Arrange taxa by column values
#'
#' @description
#' Reorder taxa based on tax_table columns or computed values. Supports the `.`
#' pronoun to refer to the phyloseq object for sorting by abundance.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param clean_phyloseq_object if TRUE (default), the resulting phyloseq object
#' is cleaned using `clean_pq()` to remove empty taxa/samples.
#' @param ... <data-masking> Variables to sort by. Use `desc()` for descending
#'   order. Use `.` to refer to the phyloseq object.
#'
#' @return A phyloseq object with reordered taxa.
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' library(MiscMetabar)
#' # Arrange by taxonomy
#' arrange_taxa_pq(data_fungi, Phylum, Class)@tax_table[, "Phylum"]
#'
#' # Arrange by total abundance (descending)
#' arrange_taxa_pq(data_fungi, dplyr::desc(taxa_sums(.)))
#' 
#' # order of columns matters
#' dfm_arr <- arrange_taxa_pq(data_fungi_mini, Class, Genus)@tax_table[, c("Class", "Genus")]
#' arrange_taxa_pq(data_fungi_mini, Genus, Class)@tax_table[, c("Class", "Genus")]
#' 
#' 
#' 
arrange_taxa_pq <- function(physeq, ..., clean_phyloseq_object = TRUE) {
  MiscMetabar::verify_pq(physeq)

  new_physeq <- physeq

  #save the phylogenetic tree if present and reorder it later
  # set phy_tree to NULL to avoid issues during reordering
 if (!is.null(phyloseq::phy_tree(physeq, errorIfNULL = FALSE))) {
    new_phy_tree <- physeq@phy_tree 
   
   new_physeq@phy_tree <- NULL
 }

  tar <- taxa_are_rows(physeq)

  tax_df <- as.data.frame(phyloseq::tax_table(physeq))
  mask <- build_taxa_data_mask(physeq)

  quos <- rlang::enquos(...)
  # Evaluate each expression in the mask context and add to data frame
  for (i in seq_along(quos)) {
    col_name <- paste0(".sort_", i)
    tax_df[[col_name]] <- rlang::eval_tidy(quos[[i]], data = mask)
  }

  # Sort by the temporary columns
  sort_cols <- paste0(".sort_", seq_along(quos))
  tax_df <- tax_df[do.call(order, tax_df[sort_cols]), , drop = FALSE]

  # Remove temporary columns
  tax_df <- tax_df[, !names(tax_df) %in% sort_cols, drop = FALSE]

  
  new_physeq@tax_table <- phyloseq::tax_table(as.matrix(tax_df))
  if (tar) {
    new_physeq@otu_table <- phyloseq::otu_table(
      as(phyloseq::otu_table(physeq), "matrix")[rownames(tax_df), , drop = FALSE],
      taxa_are_rows = TRUE
    )
  } else {
    new_physeq@otu_table <- phyloseq::otu_table(
      as(phyloseq::otu_table(physeq), "matrix")[, rownames(tax_df), drop = FALSE],
      taxa_are_rows = FALSE
    )
  }

  if (!is.null(phyloseq::phy_tree(physeq, errorIfNULL = FALSE))) {
    new_order <- rownames(tax_df)
    # Reorder tip labels and update edge matrix indices
    old_tips <- new_phy_tree$tip.label
    n_tips <- length(old_tips)
    # Create mapping: new_index[old_index] = position in new_order
    new_indices <- match(old_tips, new_order)
    # Update edge matrix: tips are numbered 1:n_tips
    is_tip <- new_phy_tree$edge <= n_tips
    new_phy_tree$edge[is_tip] <- new_indices[new_phy_tree$edge[is_tip]]
    # Reorder tip.label to match new_order
    new_phy_tree$tip.label <- new_order
    new_physeq@phy_tree <- new_phy_tree
  }

  MiscMetabar::verify_pq(new_physeq)

  if (clean_phyloseq_object) {
    return(MiscMetabar::clean_pq(new_physeq, silent = TRUE))
  } else {
    return(new_physeq)
  }
}

################################################################################
#' Rename columns in tax_table
#'
#' @description
#' Rename columns (taxonomic ranks) in tax_table using tidyselect semantics.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param ... <tidy-select> Name-value pairs where the name is the new name
#'   and the value is the old name. Use `new_name = old_name` syntax.
#'
#' @return A phyloseq object with renamed tax_table columns.
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' library(MiscMetabar)
#' # Rename a single rank
#' rename_taxa_pq(data_fungi, tax_kingdom = Kingdom)
#'
#' # Rename multiple ranks
#' rename_taxa_pq(data_fungi, tax_kingdom = Kingdom, tax_phylum = Phylum)
rename_taxa_pq <- function(physeq, ...) {
  MiscMetabar::verify_pq(physeq)
  new_physeq <- physeq

  tax_df <- as.data.frame(phyloseq::tax_table(physeq))
  tax_df <- dplyr::rename(tax_df, ...)

  new_physeq@tax_table <- phyloseq::tax_table(as.matrix(tax_df))
  phyloseq::taxa_names(new_physeq@tax_table) <- phyloseq::taxa_names(physeq)

  MiscMetabar::verify_pq(new_physeq)
  new_physeq
}
