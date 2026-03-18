test_that("filter_samples_pq filters by sample_data columns", {
  # Filter by metadata
  result <- filter_samples_pq(data_fungi, Height == "Low")
  expect_s4_class(result, "phyloseq")
  expect_true(phyloseq::nsamples(result) < phyloseq::nsamples(data_fungi))
  expect_true(all(phyloseq::sample_data(result)$Height == "Low"))
})

test_that("filter_samples_pq filters using . pronoun", {
  # Filter by sample sums
  min_depth <- 5000
  result <- filter_samples_pq(data_fungi, sample_sums(.) > min_depth)
  expect_s4_class(result, "phyloseq")
  expect_true(all(phyloseq::sample_sums(result) > min_depth))
})

test_that("filter_samples_pq combines multiple conditions", {
  result <- filter_samples_pq(
    data_fungi,
    Height == "Low",
    sample_sums(.) > 5000
  )
  expect_s4_class(result, "phyloseq")
  sam <- phyloseq::sample_data(result)
  expect_true(all(sam$Height == "Low"))
  expect_true(all(phyloseq::sample_sums(result) > 5000))
})

test_that("select_samdata_pq selects one column", {
  result <- select_samdata_pq(data_fungi, Height)
  expect_s4_class(result, "phyloseq")
  expect_equal(ncol(phyloseq::sample_data(result)), 1)
  expect_equal(colnames(phyloseq::sample_data(result)), "Height")
})

test_that("select_samdata_pq selects 2 columns", {
  result <- select_samdata_pq(data_fungi, Height, Sample_names)
  expect_s4_class(result, "phyloseq")
  expect_equal(ncol(phyloseq::sample_data(result)), 2)
  expect_equal(
    colnames(phyloseq::sample_data(result)),
    c("Height", "Sample_names")
  )
})


test_that("select_samdata_pq supports tidyselect", {
  # Test range selection
  sam_cols <- colnames(phyloseq::sample_data(data_fungi))
  if (length(sam_cols) >= 2) {
    result <- select_samdata_pq(data_fungi, 1:2)
    expect_equal(ncol(phyloseq::sample_data(result)), 2)
  }
})

test_that("mutate_samdata_pq adds new columns", {
  result <- mutate_samdata_pq(data_fungi, log_depth = log(sample_sums(.)))
  expect_s4_class(result, "phyloseq")
  expect_true("log_depth" %in% colnames(phyloseq::sample_data(result)))
  expected <- unname(log(phyloseq::sample_sums(data_fungi)))
  actual <- phyloseq::sample_data(result)$log_depth
  expect_equal(actual, expected)
})

test_that("mutate_samdata_pq modifies existing columns", {
  result <- mutate_samdata_pq(data_fungi, Height = toupper(Height))
  expect_s4_class(result, "phyloseq")
  expect_true(all(
    phyloseq::sample_data(result)$Height %in% c("LOW", "MIDDLE", "HIGH", NA)
  ))
})

test_that("slice_samples_pq selects by position", {
  result <- slice_samples_pq(data_fungi, 1:5)
  expect_s4_class(result, "phyloseq")
  expect_equal(phyloseq::nsamples(result), 5)

  original_names <- phyloseq::sample_names(data_fungi)[1:5]
  expect_equal(phyloseq::sample_names(result), original_names)
})

test_that("slice_samples_pq removes by negative index", {
  n_orig <- phyloseq::nsamples(data_fungi)
  result <- slice_samples_pq(data_fungi, -(1:2))
  expect_equal(phyloseq::nsamples(result), n_orig - 2)
})

test_that("arrange_samples_pq reorders samples", {
  result <- arrange_samples_pq(data_fungi, Height)
  expect_s4_class(result, "phyloseq")
  heights <- phyloseq::sample_data(result)$Height
  expect_true(!is.unsorted(heights))
})

test_that("arrange_samples_pq works with . pronoun", {
  result <- arrange_samples_pq(data_fungi, sample_sums(.))
  expect_s4_class(result, "phyloseq")
  sums <- phyloseq::sample_sums(result)
  expect_true(!is.unsorted(sums))
})

test_that("rename_samples_pq renames columns", {
  result <- rename_samples_pq(data_fungi, sample_height = Height)
  expect_s4_class(result, "phyloseq")
  expect_true("sample_height" %in% colnames(phyloseq::sample_data(result)))
  expect_false("Height" %in% colnames(phyloseq::sample_data(result)))
})
