# Build an object with several identical "candidate" OTUs: each has one big
# count and three low counts of 1, so it is a cross-talk candidate.
# N = 1003, n = nsamp; low threshold s*N/n >> 1, so the three 1s are low.
#   fraction  fi = Xi/N        = 3/1003        = 0.0029910
#   edgar2018 fi = (Xi/N)(n/Si) = (3/1003)(n/3) = n/1003
make_est_pq <- function(ntaxa = 12, nsamp = 12) {
  mat <- matrix(
    0,
    ntaxa,
    nsamp,
    dimnames = list(paste0("otu", seq_len(ntaxa)), paste0("s", seq_len(nsamp)))
  )
  mat[, 1] <- 1000 # big "true" sample
  mat[, 2:4] <- 1 # three low cross-talk counts
  tax <- matrix("x", ntaxa, 1, dimnames = list(rownames(mat), "Kingdom"))
  sam <- data.frame(row.names = colnames(mat), d = seq_len(nsamp))
  phyloseq::phyloseq(
    phyloseq::otu_table(mat, taxa_are_rows = TRUE),
    phyloseq::tax_table(tax),
    phyloseq::sample_data(sam)
  )
}

test_that("estimate_xtalk_rate_pq matches the hand-computed median", {
  pq <- make_est_pq(ntaxa = 12, nsamp = 12)
  ff <- estimate_xtalk_rate_pq(pq, rate = "fraction", verbose = FALSE)
  fe <- estimate_xtalk_rate_pq(pq, rate = "edgar2018", verbose = FALSE)
  expect_equal(as.numeric(ff), 3 / 1003, tolerance = 1e-9)
  expect_equal(as.numeric(fe), 12 / 1003, tolerance = 1e-9)
  expect_equal(attr(ff, "n_candidates"), 12L)
})

test_that("estimate_xtalk_rate_pq returns NA when undetermined", {
  pq <- make_est_pq(ntaxa = 5, nsamp = 12) # only 5 candidates < cmin = 10
  f <- estimate_xtalk_rate_pq(pq, verbose = FALSE)
  expect_true(is.na(f))
  expect_equal(attr(f, "n_candidates"), 5L)
})

test_that("estimate_xtalk_rate_pq respects nmin (big-OTU requirement)", {
  pq <- make_est_pq(ntaxa = 12, nsamp = 12)
  # nmin above every OTU total -> no candidates -> NA
  f <- estimate_xtalk_rate_pq(pq, nmin = 1e6, verbose = FALSE)
  expect_true(is.na(f))
  expect_equal(attr(f, "n_candidates"), 0L)
})

test_that("estimate_xtalk_rate_pq is orientation-independent", {
  pq <- make_est_pq()
  pqt <- pq
  phyloseq::otu_table(pqt) <- phyloseq::otu_table(
    t(as(phyloseq::otu_table(pq), "matrix")),
    taxa_are_rows = FALSE
  )
  expect_equal(
    as.numeric(estimate_xtalk_rate_pq(pq, verbose = FALSE)),
    as.numeric(estimate_xtalk_rate_pq(pqt, verbose = FALSE))
  )
})
