test_that("plot_sample_depth_pq returns a ggplot", {
  p <- plot_sample_depth_pq(data_fungi)
  expect_s3_class(p, "ggplot")
})

test_that("plot_sample_depth_pq works with show_threshold = FALSE", {
  p <- plot_sample_depth_pq(data_fungi, show_threshold = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_sample_depth_pq works with adjusted quantiles", {
  p <- plot_sample_depth_pq(
    data_fungi,
    lower_quantile = 0.2,
    threshold_quantile = 0.1
  )
  expect_s3_class(p, "ggplot")
})

test_that("plot_sample_depth_pq errors with too few samples", {
  pq_tiny <- slice_samples_pq(data_fungi, 1:2)
  expect_error(
    plot_sample_depth_pq(pq_tiny),
    "at least 3 samples"
  )
})
