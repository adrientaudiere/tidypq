make_nc_pq <- function() {
  mutate_samdata_pq(
    data_fungi,
    is_control = sample_sums(.) < sort(sample_sums(.))[4]
  )
}

test_that("identify_contam_negcontrol_pq returns a contam_tbl with subtype", {
  pq <- make_nc_pq()
  flagged <- identify_contam_negcontrol_pq(pq, is_control, verbose = FALSE)
  expect_true(is_contam_tbl(flagged))
  expect_true("subtype" %in% names(flagged))
  if (nrow(flagged) > 0) {
    expect_true(all(flagged$method == "negcontrol"))
  }
})

test_that("identify_contam_negcontrol_pq omits sample_contaminant by default", {
  pq <- make_nc_pq()
  default <- identify_contam_negcontrol_pq(pq, is_control, verbose = FALSE)
  expect_true(all(
    default$subtype %in% c("artifact", "lab_contaminant")
  ))
  widened <- identify_contam_negcontrol_pq(
    pq,
    is_control,
    flag_categories = c("artifact", "lab_contaminant", "sample_contaminant"),
    verbose = FALSE
  )
  expect_gte(nrow(widened), nrow(default))
})

test_that("identify_contam_negcontrol_pq validates the neg_control expression", {
  pq <- make_nc_pq()
  expect_error(
    identify_contam_negcontrol_pq(pq, sample_sums(.), verbose = FALSE),
    "logical"
  )
  pq0 <- mutate_samdata_pq(data_fungi, is_control = FALSE)
  expect_error(
    identify_contam_negcontrol_pq(pq0, is_control, verbose = FALSE),
    "No samples match"
  )
})
