#' Convert a phyloseq object to a tidy tibble
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' `pq_to_tidy()` is the canonical phyloseq-to-tidy data-preparation module for
#' the pqverse. It owns the full pipeline: cleaning, flattening, optional
#' sample aggregation, optional abundance transformation, and fact/bifactor
#' resolution. Every ggplotpq plot function calls `pq_to_tidy()` instead of
#' reinventing its own data-prep pipeline.
#'
#' The function builds the tidy tibble from scratch using `tidyr::pivot_longer()`
#' (it does not wrap `phyloseq::psmelt()`). Taxonomic rank grouping happens in
#' tibble space (no `tax_glom()`). Sample aggregation happens in tibble space
#' (no `merge_samples()` round-trip).
#'
#' **Alternative transform approaches:** users may also transform the phyloseq
#' object before calling `pq_to_tidy()` via `MiscMetabar::transform_pq()`, or
#' transform the tibble after via tidyverse `mutate()` (e.g.
#' `dplyr::mutate(abundance = abundance / sum(abundance))`).
#'
#' @param physeq (phyloseq, required) A phyloseq object. Must have
#'   `sample_data`.
#' @param fact (character or NULL, default NULL) Name of a `sample_data`
#'   column used as the primary grouping factor. `NULL` defaults to `sample_id`.
#' @param bifactor (character or NULL, default NULL) Name of a `sample_data`
#'   column used as the secondary grouping factor. Must have exactly 2 levels.
#'   `NULL` defaults to `sample_id`.
#' @param merge_sample_by (character or NULL, default NULL) Name of a
#'   `sample_data` column used to aggregate (merge) samples. Abundances are
#'   summed per group; other columns take their first value. `sample_id`
#'   becomes the group label.
#' @param transform (string, function, or NULL, default NULL) Abundance
#'   transformation applied per `sample_id` after aggregation. A string is
#'   resolved via `match.fun()` (e.g. `"log1p"`, `"sqrt"`). A function is
#'   applied directly. `NULL` leaves raw counts.
#' @param ranks (character or NULL, default NULL) Taxonomic ranks to include
#'   as columns. `NULL` includes all ranks from `tax_table`. A character vector
#'   (e.g. `c("Phylum", "Genus")`) selects specific ranks. Ignored if no
#'   `tax_table`.
#' @param filter_zero (logical, default TRUE) Drop rows where `abundance == 0`.
#' @param verbose (logical, default TRUE) Print diagnostic messages.
#'
#' @return An ungrouped tibble with columns:
#'   \itemize{
#'     \item `sample_id` (character) — sample identifier or merge group label.
#'     \item `taxon_id` (character) — OTU/ASV identifier.
#'     \item `abundance` (numeric) — counts or transformed values.
#'     \item `abundance_raw` (numeric) — pre-transform counts, present only
#'       when `transform` is applied.
#'     \item `fact` (character) — resolved primary grouping column.
#'     \item `bifactor` (character) — resolved secondary grouping column.
#'     \item Rank columns (character) — one per requested rank, `NA` replaced
#'       by `"Unknown"`. Present only if `tax_table` exists.
#'     \item All other `sample_data` columns (original names and types).
#'   }
#'
#' @seealso [phyloseq::psmelt()], [MiscMetabar::psmelt_samples_pq()],
#'   [MiscMetabar::clean_pq()], [MiscMetabar::transform_pq()]
#'
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' library(MiscMetabar)
#' # Basic tidy conversion
#' tidy_df <- pq_to_tidy(data_fungi_mini)
#' head(tidy_df)
#'
#' # With fact and bifactor
#' tidy_df <- pq_to_tidy(data_fungi_mini, fact = "Height", bifactor = "Time")
#'
#' # With sample aggregation
#' tidy_df <- pq_to_tidy(data_fungi_mini, merge_sample_by = "Height")
#'
#' # With log1p transform
#' tidy_df <- pq_to_tidy(data_fungi_mini, transform = "log1p")
#'
#' # With selected ranks only
#' tidy_df <- pq_to_tidy(data_fungi_mini, ranks = c("Phylum", "Genus"))
pq_to_tidy <- function(
  physeq,
  fact = NULL,
  bifactor = NULL,
  merge_sample_by = NULL,
  transform = NULL,
  ranks = NULL,
  filter_zero = TRUE,
  verbose = TRUE
) {
  MiscMetabar::verify_pq(physeq)

  if (is.null(phyloseq::sample_data(physeq, FALSE))) {
    stop(
      "phyloseq object must have sample_data for fact/bifactor resolution. ",
      "Use phyloseq::sample_data() to add metadata."
    )
  }

  ps <- MiscMetabar::clean_pq(
    physeq,
    remove_empty_samples = TRUE,
    remove_empty_taxa = TRUE,
    silent = !verbose
  )

  if (verbose) {
    message(
      "After clean_pq: ", phyloseq::nsamples(ps), " samples, ",
      phyloseq::ntaxa(ps), " taxa."
    )
  }

  # Step 1: Extract otu_table and normalise to taxa-as-rows
  otu_mat <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) {
    otu_mat <- t(otu_mat)
  }
  otu_df <- as.data.frame(otu_mat)
  otu_df$taxon_id <- rownames(otu_df)
  rownames(otu_df) <- NULL

  # Step 2: pivot_longer to (sample_id, taxon_id, abundance)
  df <- tidyr::pivot_longer(
    otu_df,
    cols = -taxon_id,
    names_to = "sample_id",
    values_to = "abundance"
  )

  # Step 3: Join sample_data (all columns)
  sd_df <- as.data.frame(phyloseq::sample_data(ps))
  sd_df$sample_id <- rownames(sd_df)
  rownames(sd_df) <- NULL
  df <- dplyr::left_join(df, sd_df, by = "sample_id")

  # Step 4: Join tax_table (selected ranks, NA -> "Unknown")
  has_tax_table <- !is.null(phyloseq::tax_table(ps, FALSE))
  if (has_tax_table) {
    tax_df <- as.data.frame(phyloseq::tax_table(ps))
    tax_df$taxon_id <- rownames(tax_df)
    rownames(tax_df) <- NULL

    all_ranks <- setdiff(colnames(tax_df), "taxon_id")
    if (is.null(ranks)) {
      selected_ranks <- all_ranks
    } else {
      missing_ranks <- setdiff(ranks, all_ranks)
      if (length(missing_ranks) > 0) {
        stop(
          "Ranks not found in tax_table: ",
          paste(missing_ranks, collapse = ", "),
          ". Available: ", paste(all_ranks, collapse = ", ")
        )
      }
      selected_ranks <- ranks
    }

    tax_df <- tax_df[, c("taxon_id", selected_ranks), drop = FALSE]
    for (rk in selected_ranks) {
      tax_df[[rk]][is.na(tax_df[[rk]])] <- "Unknown"
      tax_df[[rk]] <- as.character(tax_df[[rk]])
    }
    df <- dplyr::left_join(df, tax_df, by = "taxon_id")
  } else if (!is.null(ranks) && verbose) {
    message("No tax_table found; ranks parameter ignored.")
  }

  # Step 5: merge_sample_by aggregation (in tibble space)
  if (!is.null(merge_sample_by)) {
    if (!merge_sample_by %in% colnames(df)) {
      stop(
        "merge_sample_by '", merge_sample_by, "' not found in sample_data. ",
        "Available: ", paste(colnames(sd_df), collapse = ", ")
      )
    }
    if (verbose) {
      message("Aggregating samples by '", merge_sample_by, "'.")
    }
    df$sample_id <- NULL
    other_cols <- setdiff(
      colnames(df),
      c(merge_sample_by, "taxon_id", "abundance")
    )
    df <- df |>
      dplyr::group_by(.data[[merge_sample_by]], taxon_id) |>
      dplyr::summarise(
        abundance = sum(abundance),
        dplyr::across(dplyr::all_of(other_cols), dplyr::first),
        .groups = "drop"
      )
    colnames(df)[colnames(df) == merge_sample_by] <- "sample_id"
  }

  # Step 6: Transform (per sample_id, after aggregation)
  if (!is.null(transform)) {
    if (is.character(transform)) {
      transform_fn <- match.fun(transform)
    } else if (is.function(transform)) {
      transform_fn <- transform
    } else {
      stop("transform must be NULL, a string (function name), or a function.")
    }
    if (verbose) {
      message("Applying transform.")
    }
    df$abundance_raw <- df$abundance
    df <- df |>
      dplyr::group_by(sample_id) |>
      dplyr::mutate(abundance = transform_fn(abundance)) |>
      dplyr::ungroup()
  }

  # Step 7: Resolve fact/bifactor
  df$fact <- if (is.null(fact)) {
    df$sample_id
  } else {
    if (!fact %in% colnames(df)) {
      stop(
        "fact '", fact, "' not found in sample_data. Available: ",
        paste(colnames(sd_df), collapse = ", ")
      )
    }
    as.character(df[[fact]])
  }

  df$bifactor <- if (is.null(bifactor)) {
    df$sample_id
  } else {
    if (!bifactor %in% colnames(df)) {
      stop(
        "bifactor '", bifactor, "' not found in sample_data. Available: ",
        paste(colnames(sd_df), collapse = ", ")
      )
    }
    bf_vals <- as.character(df[[bifactor]])
    n_levels <- length(unique(bf_vals))
    if (n_levels != 2) {
      stop(
        "bifactor '", bifactor, "' must have exactly 2 levels, found ",
        n_levels, "."
      )
    }
    bf_vals
  }

  # Step 8: Filter zero-abundance rows
  if (filter_zero) {
    n_before <- nrow(df)
    df <- dplyr::filter(df, .data$abundance > 0)
    if (verbose) {
      message("Filtered ", n_before - nrow(df), " zero-abundance rows.")
    }
  }

  # Ensure character types for sample_id and taxon_id
  df$sample_id <- as.character(df$sample_id)
  df$taxon_id <- as.character(df$taxon_id)

  if (verbose) {
    message("Returning tibble with ", nrow(df), " rows and ",
      ncol(df), " columns.")
  }

  tibble::as_tibble(df)
}
