# Tree-scale tidyverse verbs for phyloseq objects

################################################################################
#' Filter phyloseq by tree topology
#'
#' @description
#' Filter a phyloseq object to include only taxa that are present in the
#' phylogenetic tree, or prune the tree to match the taxa in the phyloseq.
#' Can also filter based on tree properties like tip labels matching a pattern.
#'
#' @param physeq (phyloseq, required) A phyloseq object with a phy_tree slot.
#' @param taxa (character, optional) Character vector of taxa names to keep.
#'   If NULL, keeps all taxa that are present in both the OTU table and tree.
#' @param pattern (character, optional) A regular expression pattern to match
#'   against tip labels. Only tips matching the pattern are kept.
#' @param invert (logical, default FALSE) If TRUE and pattern is provided,
#'   keep tips that do NOT match the pattern.
#' @param clean_phyloseq_object if TRUE (default), the resulting phyloseq object
#' is cleaned using `clean_pq()` to remove empty taxa/samples.
#' @return A phyloseq object with filtered tree and matching taxa.
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' \dontrun{
#' library(MiscMetabar)
#' # Synchronize tree with OTU table (remove taxa not in tree)
#' filter_tree_pq(physeq_with_tree)
#'
#' # Keep only specific taxa in tree
#' filter_tree_pq(physeq_with_tree, taxa = c("ASV1", "ASV2", "ASV3"))
#'
#' # Filter by tip label pattern
#' filter_tree_pq(physeq_with_tree, pattern = "^ASV")
#' }
filter_tree_pq <- function(
  physeq,
  taxa = NULL,
  pattern = NULL,
  invert = FALSE,
  clean_phyloseq_object = TRUE
) {
  MiscMetabar::verify_pq(physeq)

  if (is.null(physeq@phy_tree)) {
    stop(
      "The phyloseq object does not contain a phylogenetic tree (phy_tree slot is NULL)."
    )
  }

  tree <- phyloseq::phy_tree(physeq)
  tip_labels <- tree$tip.label

  if (!is.null(pattern)) {
    # Filter by pattern matching on tip labels
    matches <- grepl(pattern, tip_labels)
    if (invert) {
      matches <- !matches
    }
    tips_to_keep <- tip_labels[matches]
  } else if (!is.null(taxa)) {
    # Filter by explicit taxa list
    tips_to_keep <- intersect(taxa, tip_labels)
  } else {
    # Default: synchronize tree with OTU table
    tips_to_keep <- intersect(phyloseq::taxa_names(physeq), tip_labels)
  }

  if (length(tips_to_keep) == 0) {
    stop(
      "No taxa remain after filtering. Check that taxa names match tree tip labels."
    )
  }

  # Prune tree to keep only selected tips
  new_tree <- ape::keep.tip(tree, tips_to_keep)

  # Prune physeq to keep only taxa in the pruned tree
  new_physeq <- phyloseq::prune_taxa(tips_to_keep, physeq)
  new_physeq@phy_tree <- new_tree

  if (clean_phyloseq_object) {
    return(MiscMetabar::clean_pq(new_physeq, silent = TRUE))
  } else {
    return(new_physeq)
  }
}
