test_that("identify_contam_blocklist_pq flags taxa in the blocklist", {
  flagged <- identify_contam_blocklist_pq(
    data_fungi,
    extra_genera = "Mortierella",
    verbose = FALSE
  )
  expect_true(is_contam_tbl(flagged))
  expect_true(all(flagged$method == "blocklist"))
  expect_true(nrow(flagged) > 0)
  expect_true(all(tolower(flagged$genus) == "mortierella"))
  expect_true(all(c("genus", "total_reads", "prevalence") %in% names(flagged)))
})

test_that("identify_contam_blocklist_pq returns empty when nothing matches", {
  flagged <- identify_contam_blocklist_pq(data_fungi, verbose = FALSE)
  expect_true(is_contam_tbl(flagged))
  expect_equal(nrow(flagged), 0)
})

test_that("identify_contam_blocklist_pq errors on an absent rank", {
  expect_error(
    identify_contam_blocklist_pq(data_fungi, tax_rank = "NotARank"),
    "not found"
  )
})

test_that("identify_contam_blocklist_pq respects ignore_case", {
  hit <- identify_contam_blocklist_pq(
    data_fungi,
    genus_blocklist = "mortierella",
    ignore_case = TRUE,
    verbose = FALSE
  )
  miss <- identify_contam_blocklist_pq(
    data_fungi,
    genus_blocklist = "mortierella",
    ignore_case = FALSE,
    verbose = FALSE
  )
  expect_gt(nrow(hit), nrow(miss))
})
