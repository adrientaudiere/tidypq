library(MiscMetabar)
data(data_fungi)

################################################################################
# contam_corr_pq tests
################################################################################

test_that("contam_corr_pq returns a list with expected elements", {
  res <- contam_corr_pq(data_fungi, plot = FALSE, verbose = FALSE)
  expect_type(res, "list")
  expect_named(res, c("stats", "contaminants", "summary"))
})

test_that("stats is a data.frame with expected columns", {
  res <- contam_corr_pq(data_fungi, verbose = FALSE)
  expect_s3_class(res$stats, "data.frame")
  expect_named(
    res$stats,
    c("taxon", "label", "total_reads", "prevalence",
      "cor_abs", "cor_rel", "slope", "p_value", "is_contaminant")
  )
})

test_that("summary is a named numeric vector with expected names", {
  res <- contam_corr_pq(data_fungi, verbose = FALSE)
  expect_type(res$summary, "double")
  expect_named(
    res$summary,
    c("samples_analyzed", "taxa_analyzed", "n_contaminants",
      "pct_contaminants", "contam_threshold")
  )
})

test_that("is_contaminant is a logical column", {
  res <- contam_corr_pq(data_fungi, verbose = FALSE)
  expect_type(res$stats$is_contaminant, "logical")
})

test_that("plot = FALSE returns no plot", {
  res <- contam_corr_pq(data_fungi, plot = FALSE, verbose = FALSE)
  expect_null(res$plot)
})

test_that("plot = TRUE returns a patchwork object", {
  res <- contam_corr_pq(data_fungi, plot = TRUE, verbose = FALSE)
  expect_s3_class(res$plot, "patchwork")
})

test_that("verbose = FALSE suppresses messages", {
  expect_silent(contam_corr_pq(data_fungi, verbose = FALSE))
})

test_that("error when min_total_reads is too high", {
  max_reads <- max(phyloseq::sample_sums(data_fungi), na.rm = TRUE)
  expect_error(
    contam_corr_pq(data_fungi, min_total_reads = max_reads + 1, verbose = FALSE),
    "Fewer than 3 samples have"
  )
})

test_that("contam_threshold affects number of flagged contaminants", {
  # Inject a fake contaminant (constant reads across samples)
  # so the test produces detectable contaminants
  otu <- as(phyloseq::otu_table(data_fungi), "matrix")
  if (!phyloseq::taxa_are_rows(data_fungi)) {
    otu <- t(otu)
  }
  fake <- matrix(
    rep(50, nrow(otu)),
    ncol = 1,
    dimnames = list(rownames(otu), "fake_contaminant")
  )
  otu <- cbind(otu, fake)
  tax <- rbind(
    as.matrix(phyloseq::tax_table(data_fungi)),
    fake_contaminant = rep("Contaminant", ncol(phyloseq::tax_table(data_fungi)))
  )
  pq <- phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = TRUE),
    phyloseq::tax_table(tax),
    phyloseq::sample_data(data_fungi)
  )

  res_lenient <- contam_corr_pq(pq, contam_threshold = -0.1, verbose = FALSE)
  res_strict  <- contam_corr_pq(pq, contam_threshold = -0.9, verbose = FALSE)

  # A more negative threshold (farther from 0) flags fewer contaminants,
  # since is_contaminant is based on cor_rel < threshold
  expect_lte(
    res_strict$summary[["n_contaminants"]],
    res_lenient$summary[["n_contaminants"]]
  )
})
