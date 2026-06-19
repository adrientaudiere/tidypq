# Contaminant detection via reads-vs-total distribution

################################################################################
#' Detect contaminants by correlating taxon reads with total sample depth
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Identifies potential contaminant taxa by examining the correlation
#' between each taxon's read count (or relative abundance) and the total
#' number of reads per sample. The method is inspired by the GRIMER tool
#' and relies on the observation that contaminant taxa, which are
#' introduced at a roughly constant level during sample processing, show
#' a **negative** correlation between their relative abundance and total
#' sample depth: their proportion is inflated in low-biomass samples.
#'
#' Conversely, true biological taxa tend to scale with total sample depth,
#' showing a **positive** correlation between absolute read count and
#' total reads.
#'
#' The function computes, for each taxon:
#' - Pearson correlation between absolute read count and total sample
#'   reads (`cor_abs`)
#' - Pearson correlation between relative abundance and total sample
#'   reads (`cor_rel`)
#' - Slope and p-value from a linear model of relative abundance ~
#'   total reads
#'
#' Taxa are flagged as potential contaminants when `cor_rel` is below a
#' configurable threshold (default: -0.5).
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param min_total_reads (numeric, default `1000`) Samples with fewer
#'   total reads are excluded from the analysis to avoid spurious
#'   correlations driven by extreme low-depth outliers.
#' @param min_prevalence (integer, default `2`) Taxa present in fewer
#'   than this many samples are excluded from the analysis.
#' @param contam_threshold (numeric, default `-0.5`) Taxa with
#'   `cor_rel` below this value are flagged as potential contaminants.
#'   Use values closer to 0 for a more conservative (fewer flagged)
#'   approach, or more negative for a more aggressive detection.
#' @param tax_rank (character, default `"Species"`) The taxonomic rank
#'   column used for labelling taxa in the plot. If the column does not
#'   exist, taxa names are used instead.
#' @param top_n (integer, default `20`) Number of top contaminant and
#'   top non-contaminant taxa to highlight in the diagnostic plot.
#' @param plot (logical, default `TRUE`) Whether to produce the
#'   diagnostic plot.
#' @param verbose (logical, default `TRUE`) Print a summary of detected
#'   contaminants.
#'
#' @return A list with the following elements:
#'   \itemize{
#'     \item \code{stats}: A data.frame with one row per taxon, columns:
#'       \code{taxon}, \code{total_reads}, \code{prevalence},
#'       \code{cor_abs}, \code{cor_rel}, \code{slope}, \code{p_value},
#'       \code{is_contaminant}.
#'     \item \code{contaminants}: Subset of \code{stats} for taxa
#'       flagged as contaminants.
#'     \item \code{plot}: A ggplot object (if \code{plot = TRUE}) with
#'       three panels:
#'       (A) histogram of \code{cor_rel} values,
#'       (B) scatter plot of \code{cor_abs} vs \code{cor_rel},
#'       (C) read-count vs total-depth for the top flagged taxa.
#'   }
#' @export
#' @author Adrien Taudière
#'
#' @seealso [neg_control_diag_pq()], [decontam_sam_control()]
#'
#' @examples
#' library(MiscMetabar)
#' data(data_fungi)
#'
#' # Inject a fake contaminant (constant reads across samples)
#' # so the example produces visible contaminant detection results
#' otu <- as(phyloseq::otu_table(data_fungi), "matrix")
#' fake <- matrix(rep(50, nrow(otu)), ncol = 1,
#'   dimnames = list(rownames(otu), "fake_contaminant"))
#' otu <- cbind(otu, fake)
#' tax <- rbind(
#'   as.matrix(phyloseq::tax_table(data_fungi)),
#'   fake_contaminant = rep("Contaminant", ncol(phyloseq::tax_table(data_fungi)))
#' )
#' data_fungi_c <- phyloseq::phyloseq(
#'   phyloseq::otu_table(otu, taxa_are_rows = FALSE),
#'   phyloseq::tax_table(tax),
#'   phyloseq::sample_data(data_fungi)
#' )
#'
#' # Detect contaminants
#' res <- contam_corr_pq(data_fungi_c)
#' res$summary
#' res$plot
#' head(res$contaminants)
#'
#' \dontrun{
#' # Adjust threshold
#'  res_strict <- contam_corr_pq(data_fungi_c, contam_threshold = -0.2, tax_rank = "Genus_species")
#'  res$plot
#' }
contam_corr_pq <- function(
  physeq,
  min_total_reads = 1000,
  min_prevalence = 2,
  contam_threshold = -0.5,
  tax_rank = "Species",
  top_n = 20,
  plot = TRUE,
  verbose = TRUE
) {
  MiscMetabar::verify_pq(physeq)

  # Get OTU matrix (taxa as rows)
  tar <- phyloseq::taxa_are_rows(physeq)
  otu <- as(phyloseq::otu_table(physeq), "matrix")
  if (!tar) {
    otu <- t(otu)
  }
  storage.mode(otu) <- "double"

  # Filter samples by minimum total reads
  sample_totals <- colSums(otu)
  keep_samples <- sample_totals >= min_total_reads
  if (sum(keep_samples) < 3) {
    stop(
      "Fewer than 3 samples have >= ",
      min_total_reads,
      " total reads. Lower min_total_reads or check your data.",
      call. = FALSE
    )
  }
  otu <- otu[, keep_samples, drop = FALSE]
  sample_totals <- colSums(otu)

  # Filter taxa by minimum prevalence
  prevalence <- rowSums(otu > 0)
  keep_taxa <- prevalence >= min_prevalence
  if (sum(keep_taxa) == 0) {
    stop(
      "No taxa have prevalence >= ",
      min_prevalence,
      ". Lower min_prevalence or check your data.",
      call. = FALSE
    )
  }
  otu <- otu[keep_taxa, , drop = FALSE]
  prevalence <- rowSums(otu > 0)

  # Compute relative abundances
  rel_abu <- sweep(otu, 2, sample_totals, "/")

  # Compute correlations for each taxon
  n_taxa <- nrow(otu)
  cor_abs <- numeric(n_taxa)
  cor_rel <- numeric(n_taxa)
  slopes <- numeric(n_taxa)
  p_values <- numeric(n_taxa)

  for (i in seq_len(n_taxa)) {
    abs_reads <- otu[i, ]
    rel_reads <- rel_abu[i, ]

    # Pearson correlation: absolute reads vs total reads
    cor_abs[i] <- stats::cor(abs_reads, sample_totals, use = "complete.obs")

    # Pearson correlation: relative abundance vs total reads
    cor_rel[i] <- stats::cor(rel_reads, sample_totals, use = "complete.obs")

    # Linear model: relative abundance ~ total reads
    fit <- tryCatch(
      stats::lm(rel_reads ~ sample_totals),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      slopes[i] <- stats::coef(fit)[2]
      p_val <- summary(fit)$coefficients[2, 4]
      p_values[i] <- p_val
    } else {
      slopes[i] <- NA_real_
      p_values[i] <- NA_real_
    }
  }

  # Build results data.frame
  taxa_labels <- rownames(otu)
  if (tax_rank %in% colnames(phyloseq::tax_table(physeq))) {
    tax_tab <- as.data.frame(phyloseq::tax_table(physeq))
    taxa_labels <- tax_tab[rownames(otu), tax_rank]
    taxa_labels[is.na(taxa_labels)] <- rownames(otu)[is.na(taxa_labels)]
  }

  stats_df <- data.frame(
    taxon = rownames(otu),
    label = taxa_labels,
    total_reads = rowSums(otu),
    prevalence = prevalence,
    cor_abs = cor_abs,
    cor_rel = cor_rel,
    slope = slopes,
    p_value = p_values,
    stringsAsFactors = FALSE
  )

  stats_df$is_contaminant <- stats_df$cor_rel < contam_threshold

  # Sort by cor_rel (most negative first = strongest contaminant signal)
  stats_df <- stats_df[order(stats_df$cor_rel), ]
  rownames(stats_df) <- NULL

  contaminants <- stats_df[stats_df$is_contaminant, , drop = FALSE]

  n_contam <- sum(stats_df$is_contaminant)
  n_total <- nrow(stats_df)

  if (verbose) {
    pct <- round(100 * n_contam / n_total, 1)
    cli::cli_bullets(c(
      "v" = "Contaminant detection summary:",
      "*" = "Samples analyzed: {.val {sum(keep_samples)}} (of {.val {ncol(otu) + sum(!keep_samples)}})",
      "*" = "Taxa analyzed: {.val {n_total}}",
      "*" = "Potential contaminants (cor_rel < {.val {contam_threshold}}): {.val {n_contam}} ({.val {pct}}%)",
      "*" = "Median cor_rel for contaminants: {.val {round(median(contaminants$cor_rel, na.rm = TRUE), 3)}}",
      "*" = "Median cor_rel for non-contaminants: {.val {round(median(stats_df$cor_rel[!stats_df$is_contaminant], na.rm = TRUE), 3)}}"
    ))
  }

  result <- list(
    stats = stats_df,
    contaminants = contaminants,
    summary = c(
      samples_analyzed = sum(keep_samples),
      taxa_analyzed = n_total,
      n_contaminants = n_contam,
      pct_contaminants = round(100 * n_contam / n_total, 1),
      contam_threshold = contam_threshold
    )
  )

  # Build diagnostic plot
  if (plot) {
    # Panel A: Histogram of cor_rel
    p_a <- ggplot2::ggplot(stats_df, ggplot2::aes(x = .data$cor_rel)) +
      ggplot2::geom_histogram(
        ggplot2::aes(fill = .data$is_contaminant),
        bins = 50,
        alpha = 0.75,
        boundary = 0
      ) +
      ggplot2::geom_vline(
        xintercept = contam_threshold,
        linetype = "dashed",
        colour = "red",
        linewidth = 0.8
      ) +
      ggplot2::scale_fill_manual(
        name = "Contaminant",
        values = c(`FALSE` = "#4C72B0", `TRUE` = "#D55E00"),
        labels = c(`FALSE` = "Non-contaminant", `TRUE` = "Contaminant")
      ) +
      ggplot2::labs(
        x = "Correlation (relative abundance vs total reads)",
        y = "Number of taxa",
        title = "A. Distribution of cor_rel"
      ) +
      ggplot2::annotate(
        "text",
        x = contam_threshold,
        y = Inf,
        vjust = 1.5,
        hjust = -0.1,
        label = paste("threshold =", contam_threshold),
        colour = "red",
        size = 3
      )

    # Panel B: Scatter of cor_abs vs cor_rel
    p_b <- ggplot2::ggplot(
      stats_df,
      ggplot2::aes(
        x = .data$cor_rel,
        y = .data$cor_abs
      )
    ) +
      ggplot2::geom_point(
        ggplot2::aes(colour = .data$is_contaminant),
        alpha = 0.5,
        size = 1.5
      ) +
      ggplot2::geom_vline(
        xintercept = contam_threshold,
        linetype = "dashed",
        colour = "red",
        linewidth = 0.6
      ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dotted",
        colour = "grey50"
      ) +
      ggplot2::scale_colour_manual(
        name = "Contaminant",
        values = c(`FALSE` = "#4C72B0", `TRUE` = "#D55E00"),
        labels = c(`FALSE` = "Non-contaminant", `TRUE` = "Contaminant")
      ) +
      ggplot2::labs(
        x = "Correlation (relative abundance vs total reads)",
        y = "Correlation (absolute reads vs total reads)",
        title = "B. cor_abs vs cor_rel"
      )

    # Panel C: Read-count vs total-depth for top flagged taxa
    top_contam <- utils::head(contaminants, top_n)
    if (nrow(top_contam) > 0) {
      top_taxa_ids <- top_contam$taxon
      plot_data <- data.frame(
        total_reads = rep(sample_totals, length(top_taxa_ids)),
        taxon_reads = as.vector(otu[top_taxa_ids, , drop = FALSE]),
        taxon = rep(top_contam$label, each = length(sample_totals)),
        cor_rel = rep(
          round(top_contam$cor_rel, 2),
          each = length(sample_totals)
        ),
        stringsAsFactors = FALSE
      )
      # Truncate long labels
      plot_data$taxon <- ifelse(
        nchar(plot_data$taxon) > 30,
        paste0(substr(plot_data$taxon, 1, 27), "..."),
        plot_data$taxon
      )
      plot_data$taxon <- factor(
        plot_data$taxon,
        levels = unique(plot_data$taxon)
      )

      p_c <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = .data$total_reads, y = .data$taxon_reads)
      ) +
        ggplot2::geom_point(alpha = 0.4, size = 1, colour = "#D55E00") +
        ggplot2::geom_smooth(
          method = "lm",
          se = FALSE,
          colour = "grey30",
          linewidth = 0.6,
          formula = y ~ x
        ) +
        ggplot2::facet_wrap(~taxon, scales = "free_y", ncol = 4) +
        ggplot2::labs(
          x = "Total reads per sample",
          y = "Taxon reads",
          title = "C. Top flagged taxa (reads vs sample depth)"
        )
    } else {
      p_c <- ggplot2::ggplot() +
        ggplot2::annotate(
          "text",
          x = 0,
          y = 0,
          label = "No contaminants flagged at current threshold."
        ) +
        ggplot2::theme_void() +
        ggplot2::labs(title = "C. Top flagged taxa")
    }

    result$plot <- (p_a | p_b) /
      p_c +
      patchwork::plot_layout(heights = c(1, 1.2)) +
      patchwork::plot_annotation(
        title = "Contaminant detection via reads-vs-total distribution",
        subtitle = sprintf(
          "%d taxa analyzed | %d flagged (%.1f%%) | threshold = %.2f",
          nrow(stats_df),
          sum(stats_df$is_contaminant),
          round(100 * sum(stats_df$is_contaminant) / nrow(stats_df), 1),
          contam_threshold
        )
      ) &
      ggplot2::theme(legend.position = "bottom")
  }

  return(result)
}
