# contam_tbl: canonical object returned by tidypq contaminant detectors

################################################################################
#' The `contam_tbl` class
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' A `contam_tbl` is the canonical object returned by every tidypq contaminant
#' detector (`identify_contam_*_pq()`). It is a tibble subclass with one row per
#' flagged taxon and at least two columns:
#'   \describe{
#'     \item{`taxon`}{Taxon name, matching `phyloseq::taxa_names()`.}
#'     \item{`method`}{Short token naming the detector that flagged the taxon
#'       (e.g. `"corr"`, `"blocklist"`, `"primer"`).}
#'   }
#' Detectors add their own method-specific evidence columns (correlations,
#' p-values, matched primers, …). A detector that flags nothing returns a
#' zero-row `contam_tbl`, never `NULL`, so outputs always compose.
#'
#' Outputs of several detectors can be row-bound with [rbind()]; evidence
#' columns absent from one detector are filled with `NA`. The combined table is
#' consumed by [filter_contam_pq()] (removal) and by [plot.contam_tbl()]
#' (diagnostics).
#'
#' @param x (data.frame, required) A data frame carrying at least `taxon` and
#'   `method` columns.
#'
#' @return `new_contam_tbl()` returns a `contam_tbl`. `is_contam_tbl()` returns
#'   a length-one logical.
#'
#' @name contam_tbl
#' @keywords internal
#' @author Adrien Taudière
NULL

#' @rdname contam_tbl
#' @export
new_contam_tbl <- function(
  x = tibble::tibble(
    taxon = character(),
    method = character()
  )
) {
  if (!is.data.frame(x)) {
    stop("`x` must be a data.frame.", call. = FALSE)
  }
  x <- tibble::as_tibble(x)
  if (!all(c("taxon", "method") %in% names(x))) {
    stop(
      "A contam_tbl must have `taxon` and `method` columns.",
      call. = FALSE
    )
  }
  if (nrow(x) > 0 && !is.character(x$taxon)) {
    x$taxon <- as.character(x$taxon)
  }
  if (nrow(x) > 0 && !is.character(x$method)) {
    x$method <- as.character(x$method)
  }
  # Put taxon and method first, evidence columns after.
  evidence <- setdiff(names(x), c("taxon", "method"))
  x <- x[, c("taxon", "method", evidence), drop = FALSE]
  tibble::new_tibble(x, nrow = nrow(x), class = "contam_tbl")
}

#' @rdname contam_tbl
#' @export
is_contam_tbl <- function(x) {
  inherits(x, "contam_tbl")
}

################################################################################
#' Print a `contam_tbl`
#'
#' @param x (contam_tbl, required) The object to print.
#' @param ... Passed to the tibble print method.
#'
#' @return `x`, invisibly.
#' @exportS3Method base::print
#' @author Adrien Taudière
print.contam_tbl <- function(x, ...) {
  n <- nrow(x)
  methods <- unique(x$method)
  if (n == 0) {
    cli::cli_alert_success("contam_tbl: no flagged taxa.")
  } else {
    cli::cli_text(
      "contam_tbl: {n} flagged taxon{?s} from method{?s} {.val {methods}}"
    )
  }
  NextMethod()
  invisible(x)
}

################################################################################
#' Combine several `contam_tbl` objects
#'
#' Row-binds detector outputs into a single `contam_tbl`. Evidence columns
#' missing from one detector are filled with `NA`.
#'
#' @param ... (contam_tbl) One or more `contam_tbl` objects.
#' @param deparse.level Ignored; present for `rbind()` compatibility.
#'
#' @return A `contam_tbl` holding the union of the flagged taxa.
#' @exportS3Method base::rbind
#' @author Adrien Taudière
rbind.contam_tbl <- function(..., deparse.level = 1) {
  parts <- list(...)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  combined <- dplyr::bind_rows(lapply(parts, tibble::as_tibble))
  new_contam_tbl(combined)
}

################################################################################
#' Diagnostic plot for a `contam_tbl`
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' A lightweight overview of which detectors flagged how many taxa. For
#' publication-grade contamination figures, richer plotting lives in the
#' `ggplotpq` package and consumes the same `contam_tbl`.
#'
#' @param object,x (contam_tbl, required) A detector output.
#' @param ... Currently ignored.
#'
#' @return A `ggplot` object (`autoplot`) or `x` invisibly (`plot`).
#' @importFrom ggplot2 autoplot
#' @exportS3Method ggplot2::autoplot
#' @author Adrien Taudière
autoplot.contam_tbl <- function(object, ...) {
  counts <- as.data.frame(
    table(method = object$method),
    stringsAsFactors = FALSE
  )
  if (nrow(counts) == 0) {
    counts <- data.frame(method = "none", Freq = 0)
  }
  ggplot2::ggplot(
    counts,
    ggplot2::aes(x = stats::reorder(.data$method, .data$Freq), y = .data$Freq)
  ) +
    ggplot2::geom_col(fill = "#D55E00", alpha = 0.85) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = "Detector method",
      y = "Flagged taxa",
      title = "Flagged taxa per contaminant detector"
    )
}

#' @rdname autoplot.contam_tbl
#' @exportS3Method base::plot
plot.contam_tbl <- function(x, ...) {
  print(autoplot(x, ...))
  invisible(x)
}
