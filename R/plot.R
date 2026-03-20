#' Plot sample depth differences to detect outliers
#'
#' Creates a diagnostic plot showing the log10 differences between consecutive
#' sorted sample sums. This helps identify samples with unusually low sequencing
#' depth by detecting large "jumps" in the distribution.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param lower_quantile (numeric, default: 0.1) Lower quantile threshold to
#'        exclude smallest differences when computing statistics.
#' @param threshold_quantile (numeric, default: 0.05) Quantile used to define
#'        the threshold for detecting large jumps.
#' @param show_threshold (logical, default: TRUE) Whether to color points based
#'        on the computed threshold and show the threshold rank in the title.
#'
#' @return A ggplot object showing:
#'
#'   - Points representing log10(difference) vs rank
#'   - Solid horizontal line at the mean log10(difference)
#'   - Dashed horizontal lines at the 5th and 95th percentiles
#'   - If `show_threshold = TRUE`, points colored by whether they pass the
#'     threshold
#'
#' @details
#' The function sorts samples by their total read count (sample sums), then

#' computes the difference between each consecutive pair. Large differences
#' indicate potential outlier samples with unusually low depth.
#'
#' The threshold is computed as the `threshold_quantile` of differences,
#' excluding the smallest `lower_quantile` of differences to avoid noise.
#' Samples with rank >= the first sample exceeding this threshold are
#' considered to have sufficient depth.
#'
#' @export
#'
#' @examples
#' library(MiscMetabar)
#' data(data_fungi)
#'
#' # Basic plot
#' plot_sample_depth_pq(data_fungi)
#'
#' # Adjust thresholds
#' plot_sample_depth_pq(data_fungi, lower_quantile = 0.2, threshold_quantile = 0.1)
#'
#' # Without threshold coloring
#' plot_sample_depth_pq(data_fungi, show_threshold = FALSE)

plot_sample_depth_pq <- function(
  physeq,
  lower_quantile = 0.1,
  threshold_quantile = 0.05,
  show_threshold = TRUE
) {
  MiscMetabar::verify_pq(physeq)

  if (phyloseq::nsamples(physeq) < 3) {
    stop("plot_sample_depth_pq() requires at least 3 samples.")
  }

  res_tib <-
    phyloseq::sample_sums(physeq) |>
    sort() |>
    tibble::as_tibble() |>
    dplyr::mutate(
      diff = .data$value - dplyr::lag(.data$value, default = .data$value[1]),
      rank = dplyr::row_number()
    )

  # Compute differences excluding the lower quantile and zeros (log10(0) = -Inf)
  diff_filtered <- res_tib$diff[
    res_tib$diff > stats::quantile(res_tib$diff, lower_quantile) &
      res_tib$diff > 0
  ]
  if (length(diff_filtered) == 0) {
    stop(
      "All sample depth differences are zero; cannot compute log10 statistics."
    )
  }
  mean_diff <- mean(log10(diff_filtered))
  lower_line <- stats::quantile(log10(diff_filtered), threshold_quantile)
  upper_line <- stats::quantile(log10(diff_filtered), 1 - threshold_quantile)

  p <- res_tib |>
    ggplot2::ggplot(ggplot2::aes(y = log10(.data$diff), x = .data$rank)) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = mean_diff) +
    ggplot2::geom_hline(yintercept = lower_line, linetype = 2) +
    ggplot2::geom_hline(yintercept = upper_line, linetype = 2) +
    ggplot2::labs(
      x = "Rank",
      y = "log10(difference in sample sums)"
    )

  if (show_threshold) {
    threshold <- stats::quantile(diff_filtered, threshold_quantile)
    threshold_rank <- min(which(res_tib$diff > threshold))
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(color = .data$rank >= threshold_rank)
      ) +
      ggplot2::scale_color_discrete(
        name = "Above threshold",
        labels = c("FALSE" = "No", "TRUE" = "Yes")
      ) +
      ggplot2::labs(
        title = sprintf(
          "Selected samples have rank >= %d",
          as.integer(threshold_rank)
        )
      ) +
      ggplot2::geom_vline(
        xintercept = threshold_rank,
        linetype = 4,
        color = "grey"
      ) +
      ggplot2::geom_point(
        data = res_tib[threshold_rank, ],
        ggplot2::aes(y = log10(.data$diff), x = .data$rank),
        color = "black",
        pch = 21,
        size = 4
      )
  }

  p
}
