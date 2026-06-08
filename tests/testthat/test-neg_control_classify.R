################################################################################
# neg_control_classify_pq tests
################################################################################

data(data_fungi)

mark_controls <- function() {
  mutate_samdata_pq(
    data_fungi,
    is_control = sample_sums(.) < sort(sample_sums(.))[4]
  )
}

test_that("neg_control_classify_pq returns the five expected elements", {
  pq <- mark_controls()
  res <- neg_control_classify_pq(pq, is_control, heatmap = FALSE)
  expect_named(
    res,
    c(
      "classification",
      "summary_plot",
      "distribution_plot",
      "heatmap_by_category",
      "ordination_plot"
    )
  )
  expect_s3_class(res$classification, "data.frame")
  expect_true(all(
    c("taxon", "total_reads", "category") %in% names(res$classification)
  ))
  expect_s3_class(res$summary_plot, "patchwork")
})

test_that("neg_control_classify_pq categories are the three expected levels", {
  pq <- mark_controls()
  res <- neg_control_classify_pq(pq, is_control, heatmap = FALSE)
  expect_setequal(
    levels(res$classification$category),
    c("artifact", "lab_contaminant", "sample_contaminant")
  )
})

test_that("neg_control_classify_pq skips the heatmap when heatmap = FALSE", {
  pq <- mark_controls()
  res <- neg_control_classify_pq(pq, is_control, heatmap = FALSE)
  expect_null(res$heatmap_by_category)
})

test_that("neg_control_classify_pq builds a wrapped_pheatmap when requested", {
  skip_if_not_installed("pheatmap")
  pq <- mark_controls()
  res <- neg_control_classify_pq(pq, is_control, heatmap = TRUE)
  expect_s3_class(res$heatmap_by_category, "wrapped_pheatmap")
})

test_that("neg_control_classify_pq errors when no sample matches", {
  expect_error(
    neg_control_classify_pq(
      data_fungi,
      rep(FALSE, nsamples(data_fungi)),
      heatmap = FALSE
    ),
    "No samples match"
  )
})
