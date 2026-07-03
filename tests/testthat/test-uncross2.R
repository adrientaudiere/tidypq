# Small deterministic object where otu1 dominates s1 with a 5-read leak in s2.
# With N of otu1 equal to 1005 and 4 samples, expected cross-talk z is 2.5125,
# so the leak cell scores about 0.24: flagged at tmin 0.1 but not at tmin 0.4.
make_xtalk_pq <- function() {
  mat <- matrix(
    c(
      1000,
      5,
      0,
      0, # otu1
      0,
      0,
      500,
      8 # otu2
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("otu1", "otu2"), c("s1", "s2", "s3", "s4"))
  )
  samp <- data.frame(row.names = colnames(mat), grp = c("a", "a", "b", "b"))
  tax <- matrix(
    c("k1", "k2"),
    nrow = 2,
    dimnames = list(rownames(mat), "Kingdom")
  )
  phyloseq::phyloseq(
    phyloseq::otu_table(mat, taxa_are_rows = TRUE),
    phyloseq::sample_data(samp),
    phyloseq::tax_table(tax)
  )
}

test_that("uncross2_pq returns a phyloseq and only removes reads", {
  pq <- make_xtalk_pq()
  res <- uncross2_pq(pq, verbose = FALSE)
  expect_s4_class(res, "phyloseq")
  expect_lte(sum(phyloseq::otu_table(res)), sum(phyloseq::otu_table(pq)))
})

test_that("uncross2_pq score matches the UNCROSS2 formula", {
  pq <- make_xtalk_pq()
  res <- uncross2_pq(
    pq,
    tmin = 0.1,
    return_scores = TRUE,
    clean_phyloseq_object = FALSE,
    verbose = FALSE
  )
  expect_named(
    res,
    c("old_physeq", "new_physeq", "scores", "tagjump", "tmin", "f")
  )
  expect_equal(dim(res$scores), c(2, 4))
  expect_true(all(res$scores >= 0 & res$scores <= 1))
  # z = 0.01 * 1005 / 4, x = 5
  expect_equal(
    res$scores["otu1", "s2"],
    2 / (1 + exp(5 / (0.01 * 1005 / 4))),
    tolerance = 1e-6
  )
  # empty cells score 1 (look entirely like cross-talk)
  expect_equal(res$scores["otu1", "s3"], 1)
})

test_that("tmin controls strong vs weak flagging on the borderline cell", {
  pq <- make_xtalk_pq()
  otu_val <- function(res, taxon, samp) {
    m <- as(phyloseq::otu_table(res$new_physeq), "matrix")
    if (!phyloseq::taxa_are_rows(res$new_physeq)) {
      m <- t(m)
    }
    m[taxon, samp]
  }

  weak <- uncross2_pq(
    pq,
    tmin = 0.1,
    return_scores = TRUE,
    clean_phyloseq_object = FALSE,
    verbose = FALSE
  )
  strong <- uncross2_pq(
    pq,
    tmin = 0.4,
    return_scores = TRUE,
    clean_phyloseq_object = FALSE,
    verbose = FALSE
  )

  # score ~0.24: flagged & removed at 0.1, kept at 0.4
  expect_true(weak$tagjump["otu1", "s2"])
  expect_false(strong$tagjump["otu1", "s2"])
  expect_equal(otu_val(weak, "otu1", "s2"), 0)
  expect_equal(otu_val(strong, "otu1", "s2"), 5)
  # the genuine dominant cell is never flagged
  expect_false(weak$tagjump["otu1", "s1"])
})

test_that("subtract method removes the expected background and floors at 0", {
  pq <- make_xtalk_pq()
  res <- uncross2_pq(
    pq,
    method = "subtract",
    verbose = FALSE,
    clean_phyloseq_object = FALSE
  )
  m <- as(phyloseq::otu_table(res), "matrix")
  if (!phyloseq::taxa_are_rows(res)) {
    m <- t(m)
  }
  expect_true(all(m >= 0))
  # otu1 z = 2.5125: every cell reduced by z; empty cells floored to 0
  z1 <- 0.01 * 1005 / 4
  expect_equal(m["otu1", "s1"], 1000 - z1, tolerance = 1e-6)
  expect_equal(m["otu1", "s2"], 5 - z1, tolerance = 1e-6)
  expect_equal(m["otu1", "s3"], 0)
})

test_that("uncross2_pq validates its arguments", {
  pq <- make_xtalk_pq()
  expect_error(uncross2_pq(pq, f = 0), "f")
  expect_error(uncross2_pq(pq, f = 1), "f")
  expect_error(uncross2_pq(pq, tmin = 0), "tmin")
  expect_error(uncross2_pq(pq, k = -1), "k")
})

test_that("uncross2_pq warns when empty samples inflate the denominator", {
  pq <- make_xtalk_pq()
  mat <- as(phyloseq::otu_table(pq), "matrix")
  mat[, "s4"] <- 0
  phyloseq::otu_table(pq) <- phyloseq::otu_table(mat, taxa_are_rows = TRUE)
  expect_message(uncross2_pq(pq, verbose = FALSE), "zero total reads")
})

test_that("uncross2_pq runs on data_fungi", {
  res <- uncross2_pq(data_fungi, verbose = FALSE)
  expect_s4_class(res, "phyloseq")
  expect_lte(
    sum(phyloseq::sample_sums(res)),
    sum(phyloseq::sample_sums(data_fungi))
  )
})

test_that("plot_uncross2_scores returns a ggplot for both grey settings", {
  res <- uncross2_pq(
    make_xtalk_pq(),
    return_scores = TRUE,
    clean_phyloseq_object = FALSE,
    verbose = FALSE
  )
  for (g in c(TRUE, FALSE)) {
    p <- plot_uncross2_scores(res, grey_non_significant = g)
    expect_true(inherits(p, "ggplot"))
    expect_no_error(ggplot2::ggplot_build(p))
  }
})

test_that("uncross2_pq f = 'auto' uses the de novo estimate", {
  # object with >= 10 candidate OTUs so the estimate is determinate
  mat <- matrix(
    0,
    12,
    12,
    dimnames = list(paste0("o", 1:12), paste0("s", 1:12))
  )
  mat[, 1] <- 1000
  mat[, 2:4] <- 1
  pq <- phyloseq::phyloseq(
    phyloseq::otu_table(mat, taxa_are_rows = TRUE),
    phyloseq::tax_table(matrix(
      "x",
      12,
      1,
      dimnames = list(rownames(mat), "Kingdom")
    )),
    phyloseq::sample_data(data.frame(row.names = colnames(mat), d = 1:12))
  )
  res <- uncross2_pq(
    pq,
    f = "auto",
    clean_phyloseq_object = FALSE,
    verbose = FALSE
  )
  expect_s4_class(res, "phyloseq")

  # data_fungi is undetermined under the paper's eq. (2) -> clear error
  expect_error(
    uncross2_pq(data_fungi, f = "auto", verbose = FALSE),
    "undetermined"
  )
})
