# Detector: contaminants via reads-vs-total-depth correlation (GRIMER-inspired)

################################################################################
#' Detect contaminants by correlating taxon reads with total sample depth
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Flags potential contaminant taxa from the correlation between each taxon's
#' relative abundance and total sample depth. Contaminants, introduced at a
#' roughly constant level during processing, show a **negative** correlation
#' (their proportion is inflated in low-biomass samples). Inspired by the GRIMER
#' tool.
#'
#' This is a pure detector: it flags taxa and returns a [contam_tbl] without
#' modifying `physeq`. Use [filter_contam_pq()] to remove the flagged taxa and
#' [plot.contam_tbl()] for a quick diagnostic.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param min_total_reads (numeric, default `1000`) Samples with fewer total
#'   reads are excluded to avoid spurious correlations from low-depth outliers.
#' @param min_prevalence (integer, default `2`) Taxa present in fewer than this
#'   many samples are excluded from the analysis.
#' @param contam_threshold (numeric, default `-0.5`) Taxa whose relative-
#'   abundance correlation (`cor_rel`) is below this value are flagged. Use
#'   values closer to 0 to flag more aggressively.
#' @param verbose (logical, default `TRUE`) Print a one-line detection summary.
#'
#' @return A [contam_tbl] with one row per flagged taxon, `method = "corr"`, and
#'   the evidence columns `total_reads`, `prevalence`, `cor_abs`, `cor_rel`,
#'   `slope`, and `p_value`. Empty when nothing is flagged.
#' @export
#' @author Adrien Taudière
#'
#' @seealso [filter_contam_pq()] to remove flagged taxa, [contam_tbl] for the
#'   return type.
#'
#' @examples
#' library(MiscMetabar)
#' \donttest{
#' # Inject a constant-reads fake contaminant so the example flags something
#' otu <- as(phyloseq::otu_table(data_fungi), "matrix")
#' fake <- matrix(rep(50, nrow(otu)), ncol = 1,
#'   dimnames = list(rownames(otu), "fake_contaminant"))
#' tax <- rbind(
#'   as.matrix(phyloseq::tax_table(data_fungi)),
#'   fake_contaminant = rep("Contaminant", ncol(phyloseq::tax_table(data_fungi)))
#' )
#' data_fungi_c <- phyloseq::phyloseq(
#'   phyloseq::otu_table(cbind(otu, fake), taxa_are_rows = FALSE),
#'   phyloseq::tax_table(tax),
#'   phyloseq::sample_data(data_fungi)
#' )
#'
#' flagged <- identify_contam_corr_pq(data_fungi_c)
#' flagged
#' }
identify_contam_corr_pq <- function(
  physeq,
  min_total_reads = 1000,
  min_prevalence = 2,
  contam_threshold = -0.5,
  verbose = TRUE
) {
  MiscMetabar::verify_pq(physeq)

  # OTU matrix with taxa as rows
  tar <- phyloseq::taxa_are_rows(physeq)
  otu <- as(phyloseq::otu_table(physeq), "matrix")
  if (!tar) {
    otu <- t(otu)
  }
  storage.mode(otu) <- "double"

  # Drop low-depth samples
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

  # Drop low-prevalence taxa
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

  rel_abu <- sweep(otu, 2, sample_totals, "/")

  n_taxa <- nrow(otu)
  cor_abs <- numeric(n_taxa)
  cor_rel <- numeric(n_taxa)
  slopes <- numeric(n_taxa)
  p_values <- numeric(n_taxa)

  for (i in seq_len(n_taxa)) {
    abs_reads <- otu[i, ]
    rel_reads <- rel_abu[i, ]
    cor_abs[i] <- stats::cor(abs_reads, sample_totals, use = "complete.obs")
    cor_rel[i] <- stats::cor(rel_reads, sample_totals, use = "complete.obs")
    fit <- tryCatch(
      stats::lm(rel_reads ~ sample_totals),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      slopes[i] <- stats::coef(fit)[2]
      p_values[i] <- summary(fit)$coefficients[2, 4]
    } else {
      slopes[i] <- NA_real_
      p_values[i] <- NA_real_
    }
  }

  stats_df <- tibble::tibble(
    taxon = rownames(otu),
    method = "corr",
    total_reads = rowSums(otu),
    prevalence = prevalence,
    cor_abs = cor_abs,
    cor_rel = cor_rel,
    slope = slopes,
    p_value = p_values
  )

  flagged <- stats_df[
    !is.na(stats_df$cor_rel) & stats_df$cor_rel < contam_threshold,
    ,
    drop = FALSE
  ]
  flagged <- flagged[order(flagged$cor_rel), , drop = FALSE]

  if (verbose) {
    pct <- round(100 * nrow(flagged) / nrow(stats_df), 1)
    cli::cli_alert_info(
      "corr: flagged {nrow(flagged)} of {nrow(stats_df)} taxa ({pct}%) at cor_rel < {.val {contam_threshold}}."
    )
  }

  new_contam_tbl(flagged)
}
