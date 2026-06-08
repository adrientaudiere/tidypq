# Classify negative-control taxa into contamination categories

################################################################################
#' Classify taxa found in negative controls into contamination categories
#'
#' @description
#' \lifecycle{experimental}
#'
#' Examine taxa detected in negative-control samples and assign each to one
#' of three contamination categories based on read abundance and occurrence
#' patterns:
#'
#' * **artifact**: very low total reads AND present in very few samples --
#'   likely sequencing noise or index hopping.
#' * **lab_contaminant**: predominantly found in negative controls relative
#'   to real samples -- likely introduced during library preparation or
#'   extraction.
#' * **sample_contaminant**: detected in negative controls but also
#'   widespread in real samples -- may represent genuine biological signal
#'   or index bleed.
#'
#' This complements [contam_corr_pq()] (a correlation-based method) and
#' [neg_control_diag_pq()] (a diagnostic figure) by producing an explicit
#' per-taxon category.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param neg_control <data-masking> An expression evaluated on
#'   sample_data that returns TRUE for negative-control samples (e.g.,
#'   `is_control == TRUE`, `sample_type == "blank"`). Use `.` to refer to
#'   the phyloseq object.
#' @param min_reads_artifact (numeric, default: 10) Maximum total reads
#'   (across all samples) for a taxon to be considered an artifact.
#' @param max_samples_artifact (integer, default: 2) Maximum total number
#'   of samples (negative + non-negative) in which a taxon may appear to be
#'   considered an artifact.
#' @param max_ratio_lab_contam (numeric, default: 0.2) Maximum ratio of
#'   non-negative-control sample occurrences to negative-control sample
#'   occurrences (`n_non_neg_samples / n_neg_samples`) for a taxon to be
#'   classified as a lab contaminant. Values near 0 indicate a taxon almost
#'   exclusively found in controls.
#' @param min_neg_samples_lab (integer, default: 2) Minimum number of
#'   negative-control samples a taxon must appear in to be eligible for the
#'   lab-contaminant category. Prevents rare noise from being mislabeled.
#' @param heatmap (logical, default: TRUE) If TRUE, build a
#'   \pkg{pheatmap} of negative-control abundances annotated by category.
#'   Requires the \pkg{pheatmap} package.
#'
#' @return A named list with five elements:
#'   \describe{
#'     \item{`classification`}{A `data.frame` with one row per taxon
#'       detected in at least one negative control. Columns include
#'       `taxon`, `total_reads`, `reads_in_neg`, `reads_in_samples`,
#'       `n_neg_samples`, `n_non_neg_samples`, `ratio_non_neg_to_neg`,
#'       `category`, and taxonomy columns if a `tax_table` is present.}
#'     \item{`summary_plot`}{A patchwork plot of the number of taxa and
#'       total reads per contamination category.}
#'     \item{`distribution_plot`}{A violin + jitter plot comparing the
#'       distribution of total reads across categories (log scale).}
#'     \item{`heatmap_by_category`}{A `wrapped_pheatmap` object (print it
#'       with `print()`) or `NULL` when `heatmap = FALSE`.}
#'     \item{`ordination_plot`}{A PCoA ordination (Bray-Curtis) with
#'       negative controls highlighted, or `NULL` if ordination fails.}
#'   }
#'
#' @seealso [contam_corr_pq()], [neg_control_diag_pq()],
#'   [neg_control_clean_pq()]
#'
#' @author Adrien Taudière
#' @export
#'
#' @examples
#' library(MiscMetabar)
#' data(data_fungi)
#' \donttest{
#' # Mark the three lowest-depth samples as mock negative controls
#' pq <- mutate_samdata_pq(
#'   data_fungi,
#'   is_control = sample_sums(.) < sort(sample_sums(.))[4]
#' )
#' res <- neg_control_classify_pq(pq, is_control, heatmap = FALSE)
#' table(res$classification$category)
#' res$summary_plot
#' }
neg_control_classify_pq <- function(
  physeq,
  neg_control,
  min_reads_artifact = 10,
  max_samples_artifact = 2,
  max_ratio_lab_contam = 0.2,
  min_neg_samples_lab = 2,
  heatmap = TRUE
) {
  MiscMetabar::verify_pq(physeq)

  if (heatmap && !requireNamespace("pheatmap", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg pheatmap} is required when {.code heatmap = TRUE}.",
      "i" = "Install it or call with {.code heatmap = FALSE}."
    ))
  }

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
  # otu: taxa x samples

  neg_idx <- which(nc_vec)
  non_neg_idx <- which(!nc_vec)

  otu_neg <- otu[, neg_idx, drop = FALSE]
  otu_non_neg <- otu[, non_neg_idx, drop = FALSE]

  taxa_in_neg <- rownames(otu_neg)[rowSums(otu_neg > 0) >= 1]

  if (length(taxa_in_neg) == 0) {
    cli::cli_warn("No taxa found in negative-control samples.")
    return(list(
      classification = data.frame(),
      summary_plot = ggplot2::ggplot(),
      distribution_plot = ggplot2::ggplot(),
      heatmap_by_category = NULL,
      ordination_plot = NULL
    ))
  }

  clf <- data.frame(
    taxon = taxa_in_neg,
    total_reads = rowSums(otu[taxa_in_neg, , drop = FALSE]),
    reads_in_neg = rowSums(otu_neg[taxa_in_neg, , drop = FALSE]),
    reads_in_samples = rowSums(otu_non_neg[taxa_in_neg, , drop = FALSE]),
    n_neg_samples = rowSums(otu_neg[taxa_in_neg, , drop = FALSE] > 0),
    n_non_neg_samples = rowSums(otu_non_neg[taxa_in_neg, , drop = FALSE] > 0),
    stringsAsFactors = FALSE
  )
  clf$ratio_non_neg_to_neg <- clf$n_non_neg_samples / clf$n_neg_samples

  clf$category <- NA_character_

  is_artifact <- clf$total_reads <= min_reads_artifact &
    (clf$n_neg_samples + clf$n_non_neg_samples) <= max_samples_artifact
  clf$category[is_artifact] <- "artifact"

  remaining <- is.na(clf$category)
  is_lab <- remaining &
    clf$ratio_non_neg_to_neg <= max_ratio_lab_contam &
    clf$n_neg_samples >= min_neg_samples_lab
  clf$category[is_lab] <- "lab_contaminant"

  clf$category[is.na(clf$category)] <- "sample_contaminant"
  clf$category <- factor(
    clf$category,
    levels = c("artifact", "lab_contaminant", "sample_contaminant")
  )

  tax <- tryCatch(
    as.data.frame(phyloseq::tax_table(physeq)[taxa_in_neg, ]),
    error = function(e) {
      NULL
    }
  )
  if (!is.null(tax)) {
    tax$taxon <- rownames(tax)
    clf <- dplyr::left_join(clf, tax, by = "taxon")
  }

  clf <- clf[order(clf$category, -clf$total_reads), ]
  rownames(clf) <- NULL

  # --- Summary plot: taxa count + total reads per category ---
  cat_summary <- clf |>
    dplyr::group_by(.data$category) |>
    dplyr::summarise(
      n_taxa = dplyr::n(),
      total_reads = sum(.data$total_reads),
      .groups = "drop"
    )

  p_count <- ggplot2::ggplot(
    cat_summary,
    ggplot2::aes(x = .data$category, y = .data$n_taxa, fill = .data$category)
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = .data$n_taxa), vjust = -0.3) +
    ggplot2::labs(
      title = "Number of taxa per category",
      x = NULL,
      y = "Number of taxa"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  p_reads <- ggplot2::ggplot(
    cat_summary,
    ggplot2::aes(
      x = .data$category,
      y = .data$total_reads,
      fill = .data$category
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = .data$total_reads), vjust = -0.3) +
    ggplot2::labs(
      title = "Total reads per category",
      x = NULL,
      y = "Total reads"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  summary_plot <- patchwork::wrap_plots(p_count, p_reads, ncol = 2)

  # --- Distribution plot: reads per taxon by category (+ non-contaminants) ---
  taxa_not_in_neg <- setdiff(rownames(otu), taxa_in_neg)
  cat_levels <- c(
    "artifact",
    "lab_contaminant",
    "sample_contaminant",
    "non_contaminant"
  )
  if (length(taxa_not_in_neg) > 0) {
    non_contam_df <- data.frame(
      taxon = taxa_not_in_neg,
      total_reads = rowSums(otu[taxa_not_in_neg, , drop = FALSE]),
      category = factor("non_contaminant", levels = cat_levels),
      stringsAsFactors = FALSE
    )
    distrib_df <- rbind(
      clf[, c("taxon", "total_reads", "category")],
      non_contam_df
    )
    distrib_df$category <- factor(distrib_df$category, levels = cat_levels)
  } else {
    distrib_df <- clf[, c("taxon", "total_reads", "category")]
  }

  distribution_plot <- ggplot2::ggplot(
    distrib_df,
    ggplot2::aes(
      x = .data$category,
      y = .data$total_reads + 1,
      fill = .data$category
    )
  ) +
    ggplot2::geom_violin(alpha = 0.5, scale = "width") +
    ggplot2::geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_fill_manual(
      values = c(
        artifact = "#66c2a5",
        lab_contaminant = "#fc8d62",
        sample_contaminant = "#8da0cb",
        non_contaminant = "#a6d854"
      )
    ) +
    ggplot2::labs(
      title = "Distribution of total reads per taxon by category",
      x = NULL,
      y = "Total reads (log scale)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  # --- Heatmap of log abundances in neg controls, annotated by category ---
  heatmap_wrapped <- NULL
  if (heatmap) {
    mat <- otu_neg[taxa_in_neg, , drop = FALSE]
    mat_log <- log1p(mat)
    n_taxa_per_sample <- colSums(mat > 0)
    colnames(mat_log) <- paste0(
      colnames(mat_log),
      " (",
      n_taxa_per_sample,
      ")"
    )

    annotation_row <- data.frame(
      category = clf$category,
      row.names = clf$taxon
    )
    annotation_row <- annotation_row[rownames(mat_log), , drop = FALSE]
    ann_colors <- list(
      category = c(
        artifact = "#66c2a5",
        lab_contaminant = "#fc8d62",
        sample_contaminant = "#8da0cb"
      )
    )

    heatmap_obj <- pheatmap::pheatmap(
      mat_log,
      cluster_rows = nrow(mat_log) > 1,
      cluster_cols = ncol(mat_log) > 1,
      show_rownames = nrow(mat_log) <= 60,
      show_colnames = TRUE,
      color = grDevices::colorRampPalette(c("white", "firebrick3"))(100),
      annotation_row = annotation_row,
      annotation_colors = ann_colors,
      main = "Neg. control abundances (log) by contamination category",
      silent = TRUE
    )
    heatmap_wrapped <- structure(
      list(gtable = heatmap_obj$gtable),
      class = "wrapped_pheatmap"
    )
  }

  # --- Ordination highlighting negative controls ---
  ordination_plot <- tryCatch(
    {
      physeq_clean <- phyloseq::prune_taxa(
        phyloseq::taxa_sums(physeq) > 0,
        physeq
      )
      physeq_clean <- phyloseq::prune_samples(
        phyloseq::sample_sums(physeq_clean) > 0,
        physeq_clean
      )
      physeq_log <- phyloseq::transform_sample_counts(
        physeq_clean,
        function(x) {
          log1p(x)
        }
      )
      ord <- phyloseq::ordinate(physeq_log, method = "MDS", distance = "bray")

      samp_df <- data.frame(
        phyloseq::sample_data(physeq_clean),
        check.names = FALSE
      )
      samp_df$neg_control <- nc_vec[match(
        phyloseq::sample_names(physeq_clean),
        phyloseq::sample_names(physeq)
      )]
      samp_df$sample_label <- ifelse(
        samp_df$neg_control,
        phyloseq::sample_names(physeq_clean),
        NA_character_
      )

      coords <- ord$vectors[, 1:2]
      colnames(coords) <- c("Axis1", "Axis2")
      samp_df <- cbind(samp_df, coords[phyloseq::sample_names(physeq_clean), ])

      eig <- ord$values$Eigenvalues
      var_pct <- round(100 * eig / sum(eig), 1)

      ggplot2::ggplot(samp_df, ggplot2::aes(x = .data$Axis1, y = .data$Axis2)) +
        ggplot2::geom_point(
          ggplot2::aes(
            color = .data$neg_control,
            shape = .data$neg_control
          ),
          size = 3,
          alpha = 0.8
        ) +
        ggplot2::geom_text(
          ggplot2::aes(label = .data$sample_label),
          vjust = -0.8,
          size = 3,
          na.rm = TRUE
        ) +
        ggplot2::scale_color_manual(
          values = c(`FALSE` = "grey60", `TRUE` = "red"),
          labels = c("Sample", "Negative control")
        ) +
        ggplot2::scale_shape_manual(
          values = c(`FALSE` = 16, `TRUE` = 17),
          labels = c("Sample", "Negative control")
        ) +
        ggplot2::labs(
          title = "PCoA ordination (Bray-Curtis) -- negative controls highlighted",
          x = paste0("Axis 1 (", var_pct[1], "%)"),
          y = paste0("Axis 2 (", var_pct[2], "%)"),
          color = NULL,
          shape = NULL
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "bottom")
    },
    error = function(e) {
      cli::cli_warn("Could not compute ordination plot: {e$message}")
      NULL
    }
  )

  list(
    classification = clf,
    summary_plot = summary_plot,
    distribution_plot = distribution_plot,
    heatmap_by_category = heatmap_wrapped,
    ordination_plot = ordination_plot
  )
}

################################################################################
#' Print a wrapped pheatmap object
#'
#' Draws the `gtable` stored in a `wrapped_pheatmap` object (as returned in
#' the `heatmap_by_category` slot of [neg_control_classify_pq()]).
#'
#' @param x A `wrapped_pheatmap` object.
#' @param ... Ignored, present for S3 compatibility.
#'
#' @return `x`, invisibly. Called for the side effect of drawing the heatmap.
#'
#' @author Adrien Taudière
#' @export
print.wrapped_pheatmap <- function(x, ...) {
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  invisible(x)
}
