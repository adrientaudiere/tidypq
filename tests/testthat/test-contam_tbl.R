test_that("new_contam_tbl builds a valid contam_tbl", {
  x <- new_contam_tbl(tibble::tibble(
    method = "corr",
    taxon = c("ASV1", "ASV2"),
    score = c(0.1, 0.2)
  ))
  expect_true(is_contam_tbl(x))
  expect_s3_class(x, "tbl_df")
  expect_identical(names(x)[1:2], c("taxon", "method"))
  expect_equal(nrow(x), 2)
})

test_that("new_contam_tbl defaults to an empty table", {
  e <- new_contam_tbl()
  expect_true(is_contam_tbl(e))
  expect_equal(nrow(e), 0)
  expect_true(all(c("taxon", "method") %in% names(e)))
})

test_that("new_contam_tbl requires taxon and method", {
  expect_error(
    new_contam_tbl(tibble::tibble(taxon = "ASV1")),
    "taxon.*method|method"
  )
  expect_error(new_contam_tbl("not a data frame"), "data.frame")
})

test_that("is_contam_tbl is FALSE for plain objects", {
  expect_false(is_contam_tbl(data.frame(taxon = "x", method = "y")))
  expect_false(is_contam_tbl(letters))
})

test_that("rbind.contam_tbl preserves class and fills missing evidence", {
  a <- new_contam_tbl(tibble::tibble(
    taxon = "ASV1",
    method = "corr",
    cor = -0.7
  ))
  b <- new_contam_tbl(tibble::tibble(
    taxon = "ASV2",
    method = "blocklist",
    genus = "Afipia"
  ))
  combined <- rbind(a, b)
  expect_true(is_contam_tbl(combined))
  expect_equal(nrow(combined), 2)
  expect_true(all(c("cor", "genus") %in% names(combined)))
  expect_true(is.na(combined$genus[combined$method == "corr"]))
})

test_that("print.contam_tbl returns its input invisibly", {
  x <- new_contam_tbl(tibble::tibble(taxon = "ASV1", method = "corr"))
  expect_output(out <- withVisible(print(x)))
  expect_false(out$visible)
  expect_identical(out$value, x)
})

test_that("autoplot.contam_tbl returns a ggplot", {
  x <- new_contam_tbl(tibble::tibble(
    taxon = c("ASV1", "ASV2"),
    method = c("corr", "blocklist")
  ))
  expect_s3_class(ggplot2::autoplot(x), "ggplot")
})
