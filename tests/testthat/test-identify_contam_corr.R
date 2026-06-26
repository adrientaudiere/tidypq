skip_on_cran()

make_fake_contam_pq <- function() {
  otu <- as(phyloseq::otu_table(data_fungi), "matrix")
  fake <- matrix(
    rep(50, nrow(otu)),
    ncol = 1,
    dimnames = list(rownames(otu), "fake_contaminant")
  )
  tax <- rbind(
    as.matrix(phyloseq::tax_table(data_fungi)),
    fake_contaminant = rep("Contaminant", ncol(phyloseq::tax_table(data_fungi)))
  )
  phyloseq::phyloseq(
    phyloseq::otu_table(cbind(otu, fake), taxa_are_rows = FALSE),
    phyloseq::tax_table(tax),
    phyloseq::sample_data(data_fungi)
  )
}

test_that("identify_contam_corr_pq returns a contam_tbl flagging the contaminant", {
  pq <- make_fake_contam_pq()
  flagged <- identify_contam_corr_pq(pq, verbose = FALSE)
  expect_true(is_contam_tbl(flagged))
  expect_true(all(flagged$method == "corr"))
  expect_true("fake_contaminant" %in% flagged$taxon)
  expect_true(all(
    c("cor_rel", "cor_abs", "slope", "p_value") %in% names(flagged)
  ))
})

test_that("identify_contam_corr_pq threshold controls the flag count", {
  pq <- make_fake_contam_pq()
  lax <- identify_contam_corr_pq(pq, contam_threshold = -0.1, verbose = FALSE)
  strict <- identify_contam_corr_pq(
    pq,
    contam_threshold = -0.9,
    verbose = FALSE
  )
  expect_gte(nrow(lax), nrow(strict))
})

test_that("identify_contam_corr_pq errors with too few deep samples", {
  expect_error(
    identify_contam_corr_pq(
      data_fungi,
      min_total_reads = 1e12,
      verbose = FALSE
    ),
    "Fewer than 3 samples"
  )
})
