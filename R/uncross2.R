# Filter: remove UNCROSS2-scored cross-talk / tag-jumps from an OTU table

################################################################################
#' Remove cross-talk (tag-jumps) using the UNCROSS2 score
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Detect and remove per-cell cross-talk (also called tag-jumps)
#' in an OTU table using Edgar's UNCROSS2 score (Edgar 2018,
#' \doi{10.1101/400762}). Unlike tail-trimming approaches
#' (e.g. `metabaR::tagjumpslayer()`), the score scales the expected leak by each
#' OTU's **total** abundance: an abundant genuine OTU is expected to leak more
#' reads into wells where it is absent than a rare one.
#'
#' For every OTU (row) an expected per-sample cross-talk count is computed as
#' `z = f * N / n`, where `N` is the OTU's total abundance across the whole
#' object and `n` is the **total number of samples** (Edgar's spreading model:
#' cross-talk leaks into all wells, mostly the zero/low ones). Each cell then
#' receives the UNCROSS2 score
#'
#' \deqn{s = \frac{2}{1 + \exp(k \cdot x / z)}}
#'
#' where `x` is the observed count in that cell. The score tends to `1` as
#' `x -> 0` (looks like cross-talk) and to `0` as `x >> z` (genuine signal).
#' Cells with `s >= tmin` are treated as cross-talk. The default `tmin = 0.1`
#' matches Edgar's "weak" (yellow) cut-off; use `0.4` for "strong" only.
#'
#' @section Relationship to Edgar's usearch binary: This function faithfully
#' implements the *published* UNCROSS2 score (eq. 3), verified cell-for-cell
#' against an independent reference. It differs from the usearch v11
#' `otutab_xtalk` binary in two ways that users should know:
#'
#' * **Rate.** `f` is fixed (default `0.01`) unless `f = "auto"`, whereas the
#'   binary *estimates* it de novo. Use [estimate_xtalk_rate_pq()] (or
#'   `f = "auto"`) for a de novo estimate.
#' * **Aggressiveness.** The binary's `-otutabout` filtering is far more
#'   aggressive than eq. (3): empirically it zeroes essentially every count
#'   below the per-OTU mean (`N_i/n`), keeping only the few "big" source
#'   samples, and removed ~40x more cells than eq. (3) at the same rate on data.
#'   That behaviour is undocumented and beyond the paper, so this function does
#'   **not** reproduce the binary's filtered table; it is the gentler,
#'   paper-faithful filter.
#'
#' UNCROSS2 is defined for **integer count** tables. The usearch binary rejects
#' non-integer input (`Invalid count`); this function will process relative
#' abundances but the score's interpretation no longer holds, so prefer raw
#' counts.
#'
#' @param physeq (phyloseq, required) A phyloseq object with integer counts.
#' @param f (numeric or `"auto"`, default `0.01`) Expected cross-talk rate: the
#'   fraction of an OTU's reads assumed to be misassigned across the plate. Set
#'   `f = "auto"` to estimate it de novo with [estimate_xtalk_rate_pq()]; this
#'   errors if the estimate is undetermined.
#' @param tmin (numeric, default `0.1`) Minimum UNCROSS2 score for a cell to be
#'   flagged as cross-talk. `0.1` = weak, `0.4` = strong (Edgar 2018).
#' @param k (numeric, default `1`) Steepness of the logistic score, applied as a
#'   slope on the exponent (`exp(k * x / z)`). `k = 1` reproduces UNCROSS2;
#'   `k > 1` makes the transition sharper, `k < 1` gentler.
#' @param method (character, default `"cut"`) How flagged cross-talk is handled.
#'   * `"cut"`: set every flagged cell to `0` (UNCROSS2 behaviour).
#'   * `"subtract"`: subtract the expected count `z` from **all** cells of each
#'     OTU and floor negatives at `0` (background subtraction; ignores `tmin`).
#' @param return_scores (logical, default `FALSE`) If `TRUE`, return a list with
#'   the input phyloseq object, the cleaned phyloseq object
#'   and the score / flag matrices instead of just the
#'   phyloseq object (see Value).
#' @param clean_phyloseq_object (logical, default `TRUE`) Whether to clean the
#'   resulting object with [MiscMetabar::clean_pq()] to drop emptied
#'   taxa/samples.
#' @param verbose (logical, default `TRUE`) Whether to print a one-line summary.
#'
#' @return By default a phyloseq object with flagged cross-talk removed. If
#'   `return_scores = TRUE`, a list with:
#'   * `physeq`: the cleaned phyloseq object.
#'   * `scores`: taxa-by-samples matrix of UNCROSS2 scores.
#'   * `tagjump`: logical taxa-by-samples matrix of flagged cells (always
#'     relative to `tmin`, even when `method = "subtract"`).
#'   Matrices are oriented taxa-as-rows regardless of the input orientation.
#' @export
#' @author Adrien Taudière
#'
#' @references
#' Edgar, R.C. (2018) UNCROSS2: identification of cross-talk in 16S rRNA OTU
#' tables. *bioRxiv* 400762. \doi{10.1101/400762}
#'
#' @seealso [estimate_xtalk_rate_pq()] for the de novo rate;
#'   [decontam_control_samples_pq()] for blank-based decontamination.
#'
#' @examples
#' library(MiscMetabar)
#' data(data_fungi)
#' \donttest{
#' # Remove weak cross-talk (default)
#' unc_data_fungi <- uncross2_pq(data_fungi)
#'
#' # Estimate the rate de novo instead of assuming 1%
#' uncross2_pq(data_fungi, f = "auto")
#'
#' # Only strong cross-talk, sharper transition
#' uncross2_pq(data_fungi, tmin = 0.4, k = 2)
#'
#' # Inspect scores without committing to a threshold
#' res <- uncross2_pq(data_fungi, return_scores = TRUE)
#' summary(as.vector(res$scores))
#'
#' summary_plot_pq(res$old_physeq) / summary_plot_pq(res$new_physeq)
#'
uncross2_pq <- function(
  physeq,
  f = 0.01,
  tmin = 0.1,
  k = 1,
  method = c("cut", "subtract"),
  return_scores = FALSE,
  clean_phyloseq_object = TRUE,
  verbose = TRUE
) {
  method <- match.arg(method)
  if (identical(f, "auto")) {
    f <- tidypq:::estimate_xtalk_rate_pq(physeq, verbose = verbose)
    if (is.na(f)) {
      cli::cli_abort(c(
        "de novo cross-talk rate is undetermined for this object.",
        "i" = "Supply {.arg f} manually (e.g. {.code f = 0.01})."
      ))
    }
  }
  if (!is.numeric(f) || length(f) != 1 || f <= 0 || f >= 1) {
    cli::cli_abort(
      "{.arg f} must be a single number in (0, 1), or {.val auto}."
    )
  }
  if (!is.numeric(tmin) || length(tmin) != 1 || tmin <= 0 || tmin > 1) {
    cli::cli_abort("{.arg tmin} must be a single number in (0, 1].")
  }
  if (!is.numeric(k) || length(k) != 1 || k <= 0) {
    cli::cli_abort("{.arg k} must be a single positive number.")
  }

  MiscMetabar::verify_pq(physeq, check_order = FALSE)
  old_physeq <- physeq
  physeq <- canonicalize_pq_order(old_physeq)

  if (any(phyloseq::sample_sums(physeq) == 0)) {
    cli::cli_alert_warning(c(
      "Some samples have zero total reads: they still count in {.arg n} ",
      "(total samples), inflating the denominator and shrinking expected ",
      "cross-talk. Consider {.fn MiscMetabar::clean_pq} before scoring."
    ))
  }

  # OTU matrix, taxa as rows.
  tar <- phyloseq::taxa_are_rows(physeq)
  otu <- as(phyloseq::otu_table(physeq), "matrix")
  if (!tar) {
    otu <- t(otu)
  }
  storage.mode(otu) <- "double"

  n <- ncol(otu) # total number of samples
  N <- rowSums(otu) # per-OTU total abundance
  z <- f * N / n # expected per-sample cross-talk count (length = ntaxa)

  ratio <- otu / z
  score <- 2 / (1 + exp(k * ratio))
  score[is.nan(score)] <- 0
  tagjump <- score >= tmin

  otu_clean <- otu
  if (method == "cut") {
    otu_clean[tagjump] <- 0
  } else {
    otu_clean <- otu_clean - z
    otu_clean[otu_clean < 0] <- 0
  }

  reads_removed <- sum(otu) - sum(otu_clean)
  cells_flagged <- sum(tagjump)

  otu_out <- otu_clean
  if (!tar) {
    otu_out <- t(otu_out)
  }
  new_physeq <- physeq
  new_physeq@otu_table <- phyloseq::otu_table(otu_out, taxa_are_rows = tar)

  if (clean_phyloseq_object) {
    new_physeq <- MiscMetabar::clean_pq(new_physeq, silent = TRUE)
  }

  if (verbose) {
    cli::cli_alert_info(c(
      "uncross2 ({method}): removed {reads_removed} reads",
      " across {cells_flagged} flagged cell{?s} (occurences) ",
      "(f = {f}, tmin = {tmin}, k = {k})."
    ))
  }

  MiscMetabar::verify_pq(new_physeq)

  if (return_scores) {
    return(list(
      old_physeq = old_physeq,
      new_physeq = new_physeq,
      scores = score,
      tagjump = tagjump,
      tmin = tmin,
      f = f
    ))
  } else {
    return(new_physeq)
  }
}


################################################################################
#' Plot UNCROSS2 cross-talk scores against OTU abundance
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Diagnostic scatter plot for the output of [uncross2_pq()] run with
#' `return_scores = TRUE`. Each point is one non-zero occurrence (an OTU in a
#' sample), placed by the OTU's total abundance across the object (x) and the
#' read count of that occurrence (y), and coloured by its UNCROSS2 score.
#' Occurrences scoring at or below `tmin` (kept as genuine signal) are drawn in
#' grey; those at or above `tmin` (flagged as cross-talk) are coloured on the
#' score scale. The grey 1:1 line marks where an occurrence equals the OTU
#' total: genuine "source" samples sit near it, cross-talk sits far below. Both
#' axes are on a log10 scale.
#'
#' @param scores_list (list, required) The list returned by [uncross2_pq()]
#'   with `return_scores = TRUE`. Must contain `old_physeq` (the input object),
#'   the `scores` matrix, and the `tmin` threshold.
#' @param grey_non_significant (logical, default `TRUE`) Draw occurrences
#'   scoring at or below `tmin` in grey rather than on the colour scale.
#'
#' @return A [ggplot2::ggplot] object.
#' @export
#' @author Adrien Taudière
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join filter
#' @importFrom ggplot2 ggplot geom_point geom_abline labs
#' @importFrom ggplot2 scale_color_viridis_c scale_y_log10 scale_x_log10
#'
#' @seealso [uncross2_pq()], which produces the scored list.
#'
#' @examples
#' library(MiscMetabar)
#' data(data_fungi)
#' \donttest{
#' res <- uncross2_pq(data_fungi, return_scores = TRUE)
#' plot_uncross2_scores(res)
#' }
plot_uncross2_scores <- function(scores_list, grey_non_significant = TRUE) {
  df_otutab <- scores_list$old_physeq |>
    taxa_as_rows() |>
    otu_table() |>
    as.data.frame() |>
    tibble::rownames_to_column("taxa") |>
    pivot_longer(
      -tidyselect::all_of("taxa"),
      names_to = "sample",
      values_to = "count"
    )

  df_scores <- as.data.frame(scores_list$scores) |>
    tibble::rownames_to_column("taxa") |>
    pivot_longer(
      -tidyselect::all_of("taxa"),
      names_to = "sample",
      values_to = "scores"
    )

  ts <- taxa_sums(scores_list$old_physeq)
  df_total_sums <- tibble::tibble(
    taxa = names(ts),
    total_sums = as.numeric(ts)
  )

  df <- df_scores |>
    left_join(df_total_sums, by = "taxa") |>
    left_join(df_otutab, by = c("taxa", "sample")) |>
    filter(.data$count > 0)

  p <- ggplot(
    df,
    aes(x = .data$total_sums, y = .data$count, color = .data$scores)
  )

  if (grey_non_significant) {
    # non-significant occurrences greyed, flagged ones on the colour scale
    p <- p +
      geom_point(
        data = \(d) filter(d, .data$scores <= scores_list$tmin),
        alpha = 0.2,
        color = "grey50"
      ) +
      geom_point(data = \(d) filter(d, .data$scores >= scores_list$tmin))
  } else {
    # every occurrence on the colour scale
    p <- p + geom_point()
  }

  p +
    geom_abline(slope = 1, intercept = 0, color = "grey50") +
    labs(
      x = "Total number of sequences",
      y = "Number of sequences per samples"
    ) +
    scale_color_viridis_c(na.value = "grey50") +
    scale_y_log10() +
    scale_x_log10()
}


# De novo cross-talk rate estimator (Edgar 2018, UNCROSS2, eqs 1-2).

################################################################################
#' Estimate the cross-talk rate de novo (Edgar 2018)
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Estimate the cross-talk rate `f` used by [uncross2_pq()] directly from the
#' OTU table, without control/mock samples, following the de novo procedure of
#' Edgar (2018, \doi{10.1101/400762}).
#'
#' For each OTU `i` with total abundance `N_i`, a count `c` is deemed *low*
#' (i.e. consistent with cross-talk) when `0 < c <= s*N_i/n`, where `n` is the
#' number of samples (default `s = 0.1`, so "ten times below the mean"). Let
#' `S_i` be the number of low counts and `X_i` their sum. An OTU is a
#' **candidate** when `N_i >= nmin` and `S_i >= smin`, and its rate is
#'
#' * `rate = "edgar2018"` (paper eq. 2): `f_i = (X_i/N_i)*(n/S_i)`
#' * `rate = "fraction"`: `f_i = X_i/N_i`, the raw low-read fraction, which is
#'   what the usearch `otutab_xtalk` report prints as each OTU's rate.
#'
#' Candidates with `f_i > fmax` are rejected. If fewer than `cmin` candidates
#' remain the estimate is **undetermined** and `NA` is returned (a real UNCROSS2
#' outcome; several datasets in Edgar's Table 5 are undetermined). Otherwise the
#' returned rate is the **median** `f_i` over candidates.
#'
#' @section Agreement with usearch: This reimplementation follows the *paper*.
#' The usearch v11 binary diverges from its own description in undocumented ways
#' (an absolute low-count cap of ~8 rather than `s*N_i/n`, the `X_i/N_i` rate,
#' and an iterative candidate selection), so the returned value is of the
#' **same order of magnitude** as usearch's de novo estimate but is **not**
#' bit-exact. On integer-count data the per-OTU `X_i/N_i` reproduces usearch's
#' reported rate when the same low-count set is used; the median over the whole
#' table does not. See [uncross2_pq()] for the score-side comparison.
#'
#' @param physeq (phyloseq, required) A phyloseq object with integer counts.
#' @param s (numeric, default `0.1`) Low-count fraction of the OTU mean.
#' @param nmin (numeric, default `1000`) Minimum OTU total to be a candidate.
#' @param smin (integer, default `3`) Minimum number of low counts to be a
#'   candidate.
#' @param fmax (numeric, default `0.02`) Candidates with an estimated rate above
#'   this plausible maximum are rejected.
#' @param cmin (integer, default `10`) Minimum number of candidates required;
#'   below this the estimate is undetermined (`NA`).
#' @param rate (character, default `"edgar2018"`) Per-OTU rate formula, see
#'   Description.
#' @param verbose (logical, default `TRUE`) Print a one-line summary.
#'
#' @return A single numeric: the median de novo cross-talk rate, or `NA_real_`
#'   when undetermined. The number of candidate OTUs is attached as the
#'   `"n_candidates"` attribute.
#' @keywords internal
#' @author Adrien Taudière
#'
#' @references
#' Edgar, R.C. (2018) UNCROSS2: identification of cross-talk in 16S rRNA OTU
#' tables. *bioRxiv* 400762. \doi{10.1101/400762}
#'
#' @seealso [uncross2_pq()], which accepts `f = "auto"` to use this estimate.
#'
#' @examples
#' \dontrun{
#'   estimate_xtalk_rate_pq(data_fungi)
#' }
estimate_xtalk_rate_pq <- function(
  physeq,
  s = 0.1,
  nmin = 1000,
  smin = 3L,
  fmax = 0.02,
  cmin = 10L,
  rate = c("edgar2018", "fraction"),
  verbose = TRUE
) {
  rate <- match.arg(rate)
  MiscMetabar::verify_pq(physeq, check_order = FALSE)
  physeq <- canonicalize_pq_order(physeq)

  tar <- phyloseq::taxa_are_rows(physeq)
  otu <- as(phyloseq::otu_table(physeq), "matrix")
  if (!tar) {
    otu <- t(otu)
  }
  storage.mode(otu) <- "double"

  n <- ncol(otu)
  N <- rowSums(otu)
  thr <- s * N / n # low-count threshold per OTU (recycled down columns)
  low_mask <- otu > 0 & otu <= thr
  Si <- rowSums(low_mask)
  Xi <- rowSums(otu * low_mask)

  eligible <- N >= nmin & Si >= smin
  fi <- Xi / N # fraction (safe: eligible implies N > 0)
  if (rate == "edgar2018") {
    fi <- fi * (n / Si) # Si == 0 -> Inf, masked out below
  }
  fi[!eligible] <- NA_real_

  cand <- fi[!is.na(fi) & fi <= fmax]
  n_cand <- length(cand)

  if (n_cand < cmin) {
    if (verbose) {
      cli::cli_alert_warning(c(
        "de novo cross-talk rate undetermined: only {n_cand} candidate OTU{?s} ",
        "(< cmin = {cmin}). Supply {.arg f} manually."
      ))
    }
    return(structure(NA_real_, n_candidates = n_cand))
  }

  f <- stats::median(cand)
  if (verbose) {
    cli::cli_alert_info(
      "de novo cross-talk rate ({rate}): f = {signif(f, 4)} from {n_cand} candidate OTUs."
    )
  }
  structure(f, n_candidates = n_cand)
}
