#' Remove chimeric sequences using dada2
#'
#' Applies `dada2::removeBimeraDenovo()` to identify and remove chimeric
#' sequences from a phyloseq object based on sequence abundance patterns.
#'
#' @param physeq (phyloseq, required) A phyloseq object with a refseq slot
#'   containing DNA sequences.
#' @param method (character, default: "consensus") Method for chimera detection.
#'   Passed to `dada2::removeBimeraDenovo()`. Options: "consensus", "pooled",
#'   "per-sample".
#' @param return_a_list (logical, default: FALSE) If TRUE, returns a list with
#'   the filtered phyloseq, kept taxa names, and chimeric taxa names.
#' @param ... Additional arguments passed to `dada2::removeBimeraDenovo()`.
#'
#' @return If `return_a_list = FALSE` (default), returns a phyloseq object with
#'   chimeric sequences removed. If `return_a_list = TRUE`, returns a list with:
#'   \describe{
#'     \item{physeq}{The filtered phyloseq object}
#'     \item{kept_taxa}{Character vector of retained taxa names}
#'     \item{chimeric_taxa}{Character vector of removed chimeric taxa names}
#'   }
#'
#' @details
#' This function extracts sequences and their abundances from the phyloseq
#' object, applies dada2's de novo chimera detection algorithm, and returns
#' a pruned phyloseq object containing only non-chimeric sequences.
#'
#' The dada2 method uses sequence abundance information to identify chimeras,
#' assuming that chimeric sequences are less abundant than their parent
#' sequences.
#'
#' @export
#'
#' @seealso [MiscMetabar::chimera_removal_vs()] for vsearch-based chimera
#'   removal, [create_chimera_pq()] for creating test data with synthetic
#'   chimeras.
#'
#' @examplesIf rlang::is_installed("dada2")
#' library(MiscMetabar)
#' data(data_fungi)
#'
#' \donttest{
#' # Basic usage
#' data_nochim <- chimera_removal_dada2(data_fungi)
#'
#' # Get detailed output
#' result <- chimera_removal_dada2(data_fungi, return_a_list = TRUE)
#' cat("Removed", length(result$chimeric_taxa), "chimeric ASVs\n")
#'
#' # Use pooled method
#' data_nochim <- chimera_removal_dada2(data_fungi, method = "pooled")
#' }
chimera_removal_dada2 <- function(
  physeq,
  method = "consensus",
  return_a_list = FALSE,
  ...
) {
  MiscMetabar::verify_pq(physeq)
  if (is.null(phyloseq::refseq(physeq, errorIfNULL = FALSE))) {
    stop("phyloseq object must have a refseq slot containing DNA sequences.")
  }

  df_seq_ab <- data.frame(
    sequence = as.character(phyloseq::refseq(physeq)),
    abundance = phyloseq::taxa_sums(physeq)
  )

  # Run dada2 chimera removal
  otu_nochim <- dada2::removeBimeraDenovo(
    df_seq_ab,
    method = method,
    verbose = TRUE,
    ...
  )

  # Find which ASVs were kept (compare sequences)
  kept_seqs <- rownames(otu_nochim)
  kept_taxa <- names(phyloseq::refseq(physeq))[
    names(phyloseq::refseq(physeq)) %in% kept_seqs
  ]

  if (length(kept_taxa) == 0) {
    stop(
      "No non-chimeric ASVs detected by dada2::removeBimeraDenovo(). This is unexpected."
    )
  }

  discard_taxa <- phyloseq::taxa_names(physeq)[
    !phyloseq::taxa_names(physeq) %in% kept_taxa
  ]
  discard_taxa_collapsed <- paste(
    discard_taxa[seq_len(min(10, length(discard_taxa)))],
    collapse = ";"
  )

  message(
    "Discard ",
    phyloseq::ntaxa(physeq) - length(kept_taxa),
    " chimeric ASVs including:\n",
    discard_taxa_collapsed,
    "... \nTry return_a_list=TRUE to see all discarded ASVs in the `chimeric_taxa` element."
  )

  # Subset the phyloseq object
  physeq_nochim <- phyloseq::prune_taxa(kept_taxa, physeq)

  if (return_a_list) {
    return(list(
      physeq = physeq_nochim,
      kept_taxa = kept_taxa,
      chimeric_taxa = discard_taxa
    ))
  }
  physeq_nochim
}


#' Create a phyloseq object with synthetic chimeric sequences
#'
#' @description
#' This function creates synthetic chimeric sequences by combining parts of
#' existing sequences from a phyloseq object. Useful for benchmarking chimera
#' detection methods like [MiscMetabar::chimera_removal_vs()] or
#' [chimera_removal_dada2()].
#'
#' @param physeq (phyloseq, required) A phyloseq object with a refseq slot
#'   containing DNA sequences.
#' @param n_chimeras (integer, default: 5) Number of chimeric sequences to
#'   create.
#' @param prop_mean (numeric, default: 0.5) Mean of the normal distribution
#'   used to sample the proportion of the first parent sequence.
#'   A value of 0.5 means chimeras will be centered around 50/50 splits.
#' @param prop_sd (numeric, default: 0.15) Standard deviation of the normal
#'   distribution used to sample proportions. Higher values create more
#'   variable chimera breakpoints.
#' @param prop_min (numeric, default: 0.1) Minimum proportion threshold.
#'   Proportions below this value (or above 1 - prop_min) are resampled
#'   to ensure each parent contributes meaningfully to the chimera.
#' @param seed (integer, default: 123) Random seed for reproducibility.
#' @param median_abundance_multiplier (numeric, default: 0.1) Multiplier to set
#'   the abundance of chimeric sequences relative to the median abundance of
#'   existing sequences. A value of 0.1 means chimeras will have approximately
#'   10% of the median abundance.
#' @param min_parent_distance (numeric, default: 0.1) Minimum sequence distance
#'   (proportion of differing positions) between parent1 and parent2.
#'   If 0, chimeras can be created from very similar
#'   parents, which may be harder to detect. In some cases, with min_parent_distance = 0,
#'   you may end up with chimeras that are identical to one of the parents.
#'
#'
#' @return A list containing:
#' \describe{
#'   \item{physeq}{The new phyloseq object with added chimeric sequences}
#'   \item{chimera_names}{Character vector of chimera taxa names}
#'   \item{parent_info}{Data frame with details about each chimera:
#'     chimera name, parent1, parent2, parent_distance, prop_parent1,
#'     breakpoint, seq_length}
#'   \item{params}{List of parameters used (prop_mean, prop_sd, prop_min,
#'     min_parent_distance)}
#' }
#'
#' @export
#'
#' @seealso [MiscMetabar::chimera_removal_vs()], [chimera_removal_dada2()]
#'
#' @examplesIf rlang::is_installed("Biostrings")
#' library(MiscMetabar)
#' data(data_fungi)
#'
#' # Default: centered around 50% with some variation
#' result <- create_chimera_pq(data_fungi, n_chimeras = 40)
#' data_fungi_test <- result$physeq
#' known_chimeras <- result$chimera_names
#'
#' # View the parent information and proportions
#' print(result$parent_info)
#'
#' # More variable proportions (wider distribution)
#' result2 <- create_chimera_pq(data_fungi, n_chimeras = 40,
#'                              prop_mean = 0.5, prop_sd = 0.25)
#'
#' # Biased toward more of parent1 (e.g., 70/30 splits on average)
#' result3 <- create_chimera_pq(data_fungi, n_chimeras = 40,
#'                              prop_mean = 0.7, prop_sd = 0.1)
#'
#' # Benchmark chimera detection methods
#' if (MiscMetabar::is_vsearch_installed()) {
#'   nochim_vs <- MiscMetabar::chimera_removal_vs(data_fungi_test)
#'   detected_vs <- known_chimeras[!known_chimeras %in% phyloseq::taxa_names(nochim_vs)]
#'   cat("vsearch detected:", length(detected_vs), "/",
#'       length(known_chimeras), "chimeras\n")
#' }
#'
#' # Visualize the distribution of proportions
#' hist(result$parent_info$prop_parent1,
#'      main = "Distribution of parent1 proportions",
#'      xlab = "Proportion from parent1", xlim = c(0, 1))
#'
#' # Ensure parents are at least 15% different (more detectable chimeras)
#' result4 <- create_chimera_pq(data_fungi, n_chimeras = 40,
#'                              min_parent_distance = 0.15)
#'
#' # Disable parent distance filtering (allows similar parents)
#' result5 <- create_chimera_pq(data_fungi, n_chimeras = 40,
#'                              min_parent_distance = 0)
#'
#' @author Adrien Taudiere
create_chimera_pq <- function(
  physeq,
  n_chimeras = 5,
  prop_mean = 0.5,
  prop_sd = 0.15,
  prop_min = 0.1,
  seed = 123,
  median_abundance_multiplier = 0.1,
  min_parent_distance = 0.1
) {
  set.seed(seed)

  if (is.null(phyloseq::refseq(physeq))) {
    stop("phyloseq object must have refseq slot")
  }

  if (prop_min <= 0 || prop_min >= 0.5) {
    stop("prop_min must be between 0 and 0.5 (exclusive)")
  }

  if (min_parent_distance < 0 || min_parent_distance >= 1) {
    stop("min_parent_distance must be >= 0 and < 1")
  }

  seqs <- phyloseq::refseq(physeq)
  taxa_nms <- phyloseq::taxa_names(physeq)

  # Select parent sequences (use abundant ones)
  n_candidates <- min(20, phyloseq::ntaxa(physeq))
  abundant_idx <- order(phyloseq::taxa_sums(physeq), decreasing = TRUE)[
    seq_len(n_candidates)
  ]

  # Helper function to compute sequence distance (proportion of differing positions)
  compute_seq_distance <- function(seq1, seq2) {
    len <- min(Biostrings::width(seq1), Biostrings::width(seq2))
    s1 <- as.character(Biostrings::subseq(seq1, 1, len))
    s2 <- as.character(Biostrings::subseq(seq2, 1, len))
    chars1 <- strsplit(s1, "")[[1]]
    chars2 <- strsplit(s2, "")[[1]]
    sum(chars1 != chars2) / len
  }

  chimera_list <- vector("list", n_chimeras)
  chimera_names <- character(n_chimeras)
  parent_info_rows <- vector("list", n_chimeras)

  # Function to sample proportion with minimum threshold
  sample_proportion <- function(mean, sd, min_prop) {
    prop <- stats::rnorm(1, mean = mean, sd = sd)
    # Resample if outside valid range [min_prop, 1 - min_prop]
    max_iter <- 1000
    iter <- 0
    while ((prop < min_prop || prop > (1 - min_prop)) && iter < max_iter) {
      prop <- stats::rnorm(1, mean = mean, sd = sd)
      iter <- iter + 1
    }
    if (iter == max_iter) {
      warning("Could not sample valid proportion, using prop_mean")
      prop <- mean
    }
    prop
  }

  for (i in seq_len(n_chimeras)) {
    # Pick two different parent sequences
    p1_idx <- sample(abundant_idx, 1)
    p1 <- seqs[p1_idx]

    # Select parent2, optionally ensuring sufficient distance from parent1
    if (min_parent_distance > 0) {
      candidates <- setdiff(abundant_idx, p1_idx)
      p2_idx <- NULL
      max_attempts <- length(candidates)
      attempt <- 0

      while (is.null(p2_idx) && attempt < max_attempts) {
        candidate <- sample(candidates, 1)
        dist <- compute_seq_distance(p1, seqs[candidate])
        if (dist >= min_parent_distance) {
          p2_idx <- candidate
        } else {
          candidates <- setdiff(candidates, candidate)
        }
        attempt <- attempt + 1
      }

      if (is.null(p2_idx)) {
        warning(
          "Could not find parent2 with distance >= ",
          min_parent_distance,
          " for chimera ",
          i,
          ". Using closest available candidate."
        )
        p2_idx <- sample(setdiff(abundant_idx, p1_idx), 1)
      }
    } else {
      p2_idx <- sample(setdiff(abundant_idx, p1_idx), 1)
    }

    p2 <- seqs[p2_idx]
    parents <- c(p1_idx, p2_idx)

    # Compute final parent distance
    parent_dist <- compute_seq_distance(p1, p2)

    # Sample proportion for parent1
    prop_p1 <- sample_proportion(prop_mean, prop_sd, prop_min)

    # Calculate breakpoint based on proportion
    min_len <- min(Biostrings::width(p1), Biostrings::width(p2))
    bp <- round(min_len * prop_p1)

    # Ensure breakpoint is valid
    bp <- max(1, min(bp, min_len - 1))

    # Create chimera: prop_p1 of p1 + (1-prop_p1) of p2
    chim <- Biostrings::xscat(
      Biostrings::subseq(p1, 1, bp),
      Biostrings::subseq(p2, bp + 1, min_len)
    )

    chim_name <- paste0("CHIMERA_", i)
    chimera_list[[i]] <- chim
    chimera_names[[i]] <- chim_name

    # Store parent information
    parent_info_rows[[i]] <- data.frame(
      chimera = chim_name,
      parent1 = taxa_nms[parents[1]],
      parent2 = taxa_nms[parents[2]],
      parent_distance = round(parent_dist, 4),
      prop_parent1 = round(prop_p1, 3),
      breakpoint = bp,
      seq_length = Biostrings::width(chim),
      stringsAsFactors = FALSE
    )
  }
  parent_info <- do.call(rbind, parent_info_rows)

  # Combine chimeras
  chimera_seqs <- do.call(c, chimera_list)
  names(chimera_seqs) <- chimera_names

  # Add to refseq
  new_refseq <- c(seqs, chimera_seqs)

  # Add to OTU table with moderate abundance
  otu <- as(phyloseq::otu_table(physeq), "matrix")
  if (!phyloseq::taxa_are_rows(physeq)) {
    otu <- t(otu)
  }

  # Give chimeras ~10% of median abundance
  med_abund <- stats::median(rowSums(otu))
  chim_counts <- matrix(
    stats::rpois(
      n_chimeras * ncol(otu),
      med_abund * median_abundance_multiplier / ncol(otu)
    ),
    nrow = n_chimeras,
    dimnames = list(chimera_names, colnames(otu))
  )

  new_otu <- rbind(otu, chim_counts)

  new_taxtab <- phyloseq::tax_table(physeq)
  new_taxtab <- rbind(
    as(new_taxtab, "matrix"),
    matrix(
      NA,
      nrow = n_chimeras,
      ncol = ncol(new_taxtab),
      dimnames = list(chimera_names, colnames(new_taxtab))
    )
  )

  # Rebuild phyloseq
  new_physeq <- phyloseq::phyloseq(
    phyloseq::otu_table(new_otu, taxa_are_rows = TRUE),
    phyloseq::sample_data(physeq),
    phyloseq::tax_table(new_taxtab),
    new_refseq
  )

  # Verify chimera names are in the phyloseq
  stopifnot(all(chimera_names %in% phyloseq::taxa_names(new_physeq)))

  list(
    physeq = new_physeq,
    chimera_names = chimera_names,
    parent_info = parent_info,
    params = list(
      prop_mean = prop_mean,
      prop_sd = prop_sd,
      prop_min = prop_min,
      min_parent_distance = min_parent_distance
    )
  )
}
