test_that("filter_contam_pq removes flagged taxa by name", {
  flagged <- new_contam_tbl(tibble::tibble(
    taxon = phyloseq::taxa_names(data_fungi)[1:3],
    method = "blocklist"
  ))
  res <- filter_contam_pq(data_fungi, flagged, verbose = FALSE)
  expect_s4_class(res, "phyloseq")
  expect_false(any(flagged$taxon %in% phyloseq::taxa_names(res)))
  expect_equal(phyloseq::ntaxa(res), phyloseq::ntaxa(data_fungi) - 3)
})

test_that("filter_contam_pq is a no-op on an empty contam_tbl", {
  res <- filter_contam_pq(data_fungi, new_contam_tbl(), verbose = FALSE)
  expect_equal(phyloseq::ntaxa(res), phyloseq::ntaxa(data_fungi))
})

test_that("filter_contam_pq de-duplicates taxa flagged by several methods", {
  dup <- new_contam_tbl(tibble::tibble(
    taxon = rep(phyloseq::taxa_names(data_fungi)[1], 2),
    method = c("corr", "blocklist")
  ))
  res <- filter_contam_pq(data_fungi, dup, verbose = FALSE)
  expect_equal(phyloseq::ntaxa(res), phyloseq::ntaxa(data_fungi) - 1)
})

test_that("filter_contam_pq rejects non-contam_tbl input", {
  expect_error(
    filter_contam_pq(data_fungi, phyloseq::taxa_names(data_fungi)[1:2]),
    "contam_tbl"
  )
  expect_error(
    filter_contam_pq(data_fungi, data.frame(taxon = "x", method = "y")),
    "contam_tbl"
  )
})
