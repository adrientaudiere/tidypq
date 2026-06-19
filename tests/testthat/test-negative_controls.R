################################################################################
# neg_control_diag_pq tests
################################################################################

test_that("neg_control_diag_pq returns a patchwork object", {
  pq <- mutate_samdata_pq(
    data_fungi,
    is_control = seq_len(nsamples(.)) <= 2
  )
  fig <- suppressWarnings(neg_control_diag_pq(pq, is_control))
  expect_s3_class(fig, "patchwork")
})

test_that("neg_control_diag_pq works with log10_transform = FALSE", {
  pq <- mutate_samdata_pq(
    data_fungi,
    is_control = seq_len(nsamples(.)) <= 2
  )
  fig <- suppressWarnings(neg_control_diag_pq(pq, is_control, log10_transform = FALSE))
  expect_s3_class(fig, "patchwork")
})

test_that("neg_control_diag_pq accepts a custom palette and title", {
  pq <- mutate_samdata_pq(
    data_fungi,
    is_control = seq_len(nsamples(.)) <= 2
  )
  fig <- suppressWarnings(neg_control_diag_pq(
    pq,
    is_control,
    palette = c(`FALSE` = "#888888", `TRUE` = "#000000"),
    title = "Custom title"
  ))
  expect_s3_class(fig, "patchwork")
})

test_that("neg_control_diag_pq evaluates data-masked expressions", {
  # control_condition can use sample_data columns directly
  pq <- mutate_samdata_pq(
    data_fungi,
    is_control = seq_len(nsamples(.)) <= 2
  )
  fig <- suppressWarnings(neg_control_diag_pq(pq, is_control == TRUE))
  expect_s3_class(fig, "patchwork")
})

test_that("neg_control_diag_pq errors when no samples match", {
  pq <- mutate_samdata_pq(
    data_fungi,
    is_control = rep(FALSE, nsamples(.))
  )
  expect_error(
    neg_control_diag_pq(pq, is_control),
    "No samples match"
  )
})

test_that("neg_control_diag_pq errors when condition is not logical", {
  expect_error(
    neg_control_diag_pq(data_fungi, sample_sums(.)),
    "must evaluate to a logical vector"
  )
})

test_that("neg_control_diag_pq returns a single-panel figure when n_real == 0", {
  pq <- mutate_samdata_pq(
    data_fungi,
    is_control = rep(TRUE, nsamples(.))
  )
  suppressWarnings(expect_warning(
    fig <- neg_control_diag_pq(pq, is_control),
    "only panel A"
  ))
  expect_s3_class(fig, "patchwork")
})

################################################################################
# neg_control_clean_pq tests
################################################################################

test_that("neg_control_clean_pq returns a phyloseq object", {
  pq <- mutate_samdata_pq(
    data_fungi,
    is_control = seq_len(nsamples(.)) <= 2
  )
  result <- neg_control_clean_pq(pq, is_control, verbose = FALSE)
  expect_s4_class(result, "phyloseq")
})

test_that("neg_control_clean_pq drops NC samples by default", {
  pq <- mutate_samdata_pq(
    data_fungi,
    is_control = seq_len(nsamples(.)) <= 2
  )
  result <- neg_control_clean_pq(pq, is_control, verbose = FALSE)
  nc_in <- phyloseq::sample_names(pq)[seq_len(2)]
  expect_false(any(phyloseq::sample_names(result) %in% nc_in))
})

test_that("neg_control_clean_pq keeps NC samples when drop_neg_samples = FALSE", {
  pq <- mutate_samdata_pq(
    data_fungi,
    is_control = seq_len(nsamples(.)) <= 2
  )
  result <- neg_control_clean_pq(
    pq,
    is_control,
    drop_neg_samples = FALSE,
    verbose = FALSE
  )
  nc_in <- phyloseq::sample_names(pq)[seq_len(2)]
  expect_true(all(nc_in %in% phyloseq::sample_names(result)))
})

test_that("neg_control_clean_pq removes contaminant taxa above threshold", {
  pq <- mutate_samdata_pq(
    data_fungi,
    is_control = seq_len(nsamples(.)) <= 2
  )
  nc_only <- phyloseq::prune_samples(
    seq_len(phyloseq::nsamples(pq)) <= 2,
    pq
  )
  nc_sums <- phyloseq::taxa_sums(nc_only)
  expected_drop <- names(nc_sums)[nc_sums > 2]

  result <- neg_control_clean_pq(
    pq,
    is_control,
    min_reads = 2,
    drop_neg_samples = FALSE,
    clean_phyloseq_object = FALSE,
    verbose = FALSE
  )
  # none of the contaminants should remain
  expect_false(any(expected_drop %in% phyloseq::taxa_names(result)))
  # all non-contaminants should remain
  kept <- setdiff(phyloseq::taxa_names(pq), expected_drop)
  expect_true(all(kept %in% phyloseq::taxa_names(result)))
})

test_that("neg_control_clean_pq with min_reads = 0 drops any NC-detected taxon", {
  pq <- mutate_samdata_pq(
    data_fungi,
    is_control = seq_len(nsamples(.)) <= 2
  )
  result <- neg_control_clean_pq(
    pq,
    is_control,
    min_reads = 0,
    drop_neg_samples = FALSE,
    clean_phyloseq_object = FALSE,
    verbose = FALSE
  )
  nc_only <- phyloseq::prune_samples(
    seq_len(phyloseq::nsamples(pq)) <= 2,
    pq
  )
  detected_in_nc <- names(phyloseq::taxa_sums(nc_only))[
    phyloseq::taxa_sums(nc_only) > 0
  ]
  expect_false(any(detected_in_nc %in% phyloseq::taxa_names(result)))
})

test_that("neg_control_clean_pq with very high min_reads keeps every taxon", {
  pq <- mutate_samdata_pq(
    data_fungi,
    is_control = seq_len(nsamples(.)) <= 2
  )
  huge <- max(phyloseq::taxa_sums(pq)) + 1
  result <- neg_control_clean_pq(
    pq,
    is_control,
    min_reads = huge,
    drop_neg_samples = FALSE,
    clean_phyloseq_object = FALSE,
    verbose = FALSE
  )
  expect_equal(
    phyloseq::ntaxa(result),
    phyloseq::ntaxa(pq)
  )
})

test_that("neg_control_clean_pq emits a message when verbose = TRUE", {
  pq <- mutate_samdata_pq(
    data_fungi,
    is_control = seq_len(nsamples(.)) <= 2
  )
  expect_message(
    neg_control_clean_pq(pq, is_control, min_reads = 0, verbose = TRUE),
    "Removing"
  )
})

test_that("neg_control_clean_pq errors when no samples match", {
  pq <- mutate_samdata_pq(
    data_fungi,
    is_control = rep(FALSE, nsamples(.))
  )
  expect_error(
    neg_control_clean_pq(pq, is_control, verbose = FALSE),
    "No samples match"
  )
})

test_that("neg_control_clean_pq errors when condition is not logical", {
  expect_error(
    neg_control_clean_pq(data_fungi, sample_sums(.), verbose = FALSE),
    "must evaluate to a logical vector"
  )
})
