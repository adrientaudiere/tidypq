test_that("pq_to_tidy basic output", {
  skip_if_not_installed("MiscMetabar")
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  df <- pq_to_tidy(data_fungi_mini, verbose = FALSE)

  expect_s3_class(df, "tbl_df")
  expect_true("sample_id" %in% colnames(df))
  expect_true("taxon_id" %in% colnames(df))
  expect_true("abundance" %in% colnames(df))
  expect_true("fact" %in% colnames(df))
  expect_true("bifactor" %in% colnames(df))
  expect_false("abundance_raw" %in% colnames(df))
})

test_that("pq_to_tidy includes all sample_data columns", {
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  df <- pq_to_tidy(data_fungi_mini, verbose = FALSE)
  sd_names <- colnames(as.data.frame(
    phyloseq::sample_data(data_fungi_mini)
  ))
  expect_true(all(sd_names %in% colnames(df)))
})

test_that("pq_to_tidy includes rank columns from tax_table", {
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  df <- pq_to_tidy(data_fungi_mini, verbose = FALSE)
  tax_ranks <- colnames(phyloseq::tax_table(data_fungi_mini))
  expect_true(all(tax_ranks %in% colnames(df)))
})

test_that("pq_to_tidy psmelt parity", {
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  df <- pq_to_tidy(data_fungi_mini, verbose = FALSE, filter_zero = FALSE)
  psm <- as.data.frame(phyloseq::psmelt(data_fungi_mini))

  df_sorted <- df[order(df$sample_id, df$taxon_id), ]
  psm_sorted <- psm[order(psm$Sample, psm$OTU), ]

  expect_equal(
    sort(unique(df$sample_id)),
    sort(unique(psm_sorted$Sample))
  )
  expect_equal(
    sort(unique(df$taxon_id)),
    sort(unique(psm_sorted$OTU))
  )

  df_ab <- df_sorted$abundance
  psm_ab <- psm_sorted$Abundance
  expect_equal(df_ab, psm_ab)
})

test_that("pq_to_tidy fact defaults to sample_id", {
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  df <- pq_to_tidy(data_fungi_mini, verbose = FALSE)
  expect_equal(df$fact, df$sample_id)
})

test_that("pq_to_tidy bifactor defaults to sample_id", {
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  df <- pq_to_tidy(data_fungi_mini, verbose = FALSE)
  expect_equal(df$bifactor, df$sample_id)
})

test_that("pq_to_tidy fact resolves named column", {
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  df <- pq_to_tidy(data_fungi_mini, fact = "Height", verbose = FALSE)
  expect_equal(df$fact, as.character(df$Height))
})

test_that("pq_to_tidy bifactor enforces 2 levels", {
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  expect_error(
    pq_to_tidy(data_fungi_mini, bifactor = "Height", verbose = FALSE),
    "exactly 2 levels"
  )
})

test_that("pq_to_tidy bifactor with valid 2-level column", {
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  ps <- data_fungi_mini
  sd <- as.data.frame(phyloseq::sample_data(ps))
  sd$two_level <- factor(rep(c("A", "B"), length.out = nrow(sd)))
  phyloseq::sample_data(ps) <- phyloseq::sample_data(sd)

  df <- pq_to_tidy(ps, bifactor = "two_level", verbose = FALSE)
  expect_true(all(df$bifactor %in% c("A", "B")))
})

test_that("pq_to_tidy errors on invalid fact", {
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  expect_error(
    pq_to_tidy(data_fungi_mini, fact = "nonexistent", verbose = FALSE),
    "not found in sample_data"
  )
})

test_that("pq_to_tidy errors on missing sample_data", {
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  ps_no_sd <- phyloseq::phyloseq(
    phyloseq::otu_table(data_fungi_mini),
    phyloseq::tax_table(data_fungi_mini)
  )

  expect_error(
    pq_to_tidy(ps_no_sd, verbose = FALSE),
    "must have sample_data"
  )
})

test_that("pq_to_tidy merge_sample_by aggregation", {
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  df <- pq_to_tidy(data_fungi_mini, merge_sample_by = "Height", verbose = FALSE)

  expect_equal(
    sort(unique(df$sample_id)),
    sort(unique(as.character(
      phyloseq::sample_data(data_fungi_mini)$Height
    )))
  )

  df_agg <- dplyr::summarise(
    dplyr::group_by(df, sample_id, taxon_id),
    abundance = sum(abundance),
    .groups = "drop"
  )
  expect_equal(nrow(df), nrow(df_agg))
})

test_that("pq_to_tidy transform with function", {
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  df <- pq_to_tidy(
    data_fungi_mini,
    transform = function(x) x / sum(x),
    verbose = FALSE
  )
  expect_true("abundance_raw" %in% colnames(df))
  expect_true(all(df$abundance >= 0))
  sample_sums <- tapply(df$abundance, df$sample_id, sum)
  expect_true(all(abs(sample_sums - 1) < 1e-10))
})

test_that("pq_to_tidy transform with string via match.fun", {
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  df <- pq_to_tidy(data_fungi_mini, transform = "log1p", verbose = FALSE)
  expect_true("abundance_raw" %in% colnames(df))
  expect_true(all(df$abundance >= 0))
})

test_that("pq_to_tidy transform preserves abundance_raw", {
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  df_raw <- pq_to_tidy(data_fungi_mini, verbose = FALSE)
  df_trans <- pq_to_tidy(data_fungi_mini, transform = "log1p", verbose = FALSE)
  expect_equal(df_trans$abundance_raw, df_raw$abundance)
})

test_that("pq_to_tidy ranks filtering", {
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  all_ranks <- colnames(phyloseq::tax_table(data_fungi_mini))
  selected <- all_ranks[1:2]

  df <- pq_to_tidy(data_fungi_mini, ranks = selected, verbose = FALSE)
  expect_true(all(selected %in% colnames(df)))
  missing_ranks <- setdiff(all_ranks, selected)
  expect_false(any(missing_ranks %in% colnames(df)))
})

test_that("pq_to_tidy NA in taxonomy becomes Unknown", {
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  ps <- data_fungi_mini
  tax <- as.data.frame(phyloseq::tax_table(ps))
  tax[1, ncol(tax)] <- NA
  phyloseq::tax_table(ps) <- phyloseq::tax_table(as.matrix(tax))

  df <- pq_to_tidy(ps, verbose = FALSE)
  expect_false(any(is.na(df[[ncol(tax) + 3]])))
})

test_that("pq_to_tidy filter_zero removes zero rows", {
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  df_filtered <- pq_to_tidy(
    data_fungi_mini,
    filter_zero = TRUE,
    verbose = FALSE
  )
  df_all <- pq_to_tidy(data_fungi_mini, filter_zero = FALSE, verbose = FALSE)
  expect_true(nrow(df_filtered) < nrow(df_all))
  expect_true(all(df_filtered$abundance > 0))
})

test_that("pq_to_tidy handles taxa-as-columns", {
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  ps <- data_fungi_mini
  if (phyloseq::taxa_are_rows(ps)) {
    otu <- as(phyloseq::otu_table(ps), "matrix")
    phyloseq::otu_table(ps) <- phyloseq::otu_table(
      t(otu),
      taxa_are_rows = FALSE
    )
  }

  df <- pq_to_tidy(ps, verbose = FALSE)
  expect_s3_class(df, "tbl_df")
  expect_true(all(df$abundance >= 0))
})

test_that("pq_to_tidy returns ungrouped tibble", {
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  df <- pq_to_tidy(data_fungi_mini, verbose = FALSE)
  expect_length(dplyr::group_vars(df), 0)
})

test_that("pq_to_tidy without tax_table", {
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  ps <- data_fungi_mini
  phyloseq::tax_table(ps) <- NULL

  df <- pq_to_tidy(ps, verbose = FALSE)
  expect_s3_class(df, "tbl_df")
  tax_ranks <- colnames(phyloseq::tax_table(data_fungi_mini))
  expect_false(any(tax_ranks %in% colnames(df)))
})

test_that("pq_to_tidy verbose outputs messages", {
  data_fungi_mini <- NULL
  data("data_fungi_mini", package = "MiscMetabar", envir = environment())

  expect_message(pq_to_tidy(data_fungi_mini, verbose = TRUE))
  expect_silent(pq_to_tidy(data_fungi_mini, verbose = FALSE))
})
