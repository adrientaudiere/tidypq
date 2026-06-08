################################################################################
# contam_blocklist_pq tests
################################################################################

data(data_fungi)

test_that("known_contaminant_genera returns a non-empty character vector", {
  bl <- known_contaminant_genera()
  expect_type(bl, "character")
  expect_gt(length(bl), 50)
  expect_false(any(duplicated(bl)))
})

test_that("contam_blocklist_pq returns the expected structure", {
  res <- contam_blocklist_pq(data_fungi, verbose = FALSE)
  expect_named(res, c("stats", "contaminants", "blocklist", "summary"))
  expect_s3_class(res$stats, "data.frame")
  expect_equal(nrow(res$stats), ntaxa(data_fungi))
  expect_true(all(
    c(
      "taxon",
      "genus",
      "is_known_contaminant",
      "total_reads",
      "prevalence"
    ) %in%
      names(res$stats)
  ))
})

test_that("contam_blocklist_pq flags genera added via extra_genera", {
  genus_present <- stats::na.omit(
    as.vector(phyloseq::tax_table(data_fungi)[, "Genus"])
  )[1]
  res <- contam_blocklist_pq(
    data_fungi,
    extra_genera = genus_present,
    verbose = FALSE
  )
  expect_true(any(res$stats$is_known_contaminant))
  expect_true(genus_present %in% res$contaminants$genus)
})

test_that("contam_blocklist_pq matching is case-insensitive by default", {
  genus_present <- stats::na.omit(
    as.vector(phyloseq::tax_table(data_fungi)[, "Genus"])
  )[1]
  res_lower <- contam_blocklist_pq(
    data_fungi,
    genus_blocklist = tolower(genus_present),
    verbose = FALSE
  )
  expect_true(genus_present %in% res_lower$contaminants$genus)

  res_strict <- contam_blocklist_pq(
    data_fungi,
    genus_blocklist = tolower(genus_present),
    ignore_case = FALSE,
    verbose = FALSE
  )
  expect_false(genus_present %in% res_strict$contaminants$genus)
})

test_that("contam_blocklist_pq errors on a missing rank", {
  expect_error(
    contam_blocklist_pq(data_fungi, tax_rank = "NotARank", verbose = FALSE),
    "not found"
  )
})
