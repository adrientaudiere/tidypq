test_that("filter_occurrences_pq removes values below threshold", {
  # Remove singletons
  result <- filter_occurrences_pq(data_fungi_mini, . > 1, clean_phyloseq_object = FALSE)
  expect_s4_class(result, "phyloseq")

  otu_orig <- as(phyloseq::otu_table(data_fungi_mini), "matrix")
  otu_new <- as(phyloseq::otu_table(result), "matrix")

  # Check that singleton values are now 0
  singleton_positions <- which(otu_orig == 1, arr.ind = TRUE)
  if (nrow(singleton_positions) > 0) {
    for (i in seq_len(nrow(singleton_positions))) {
      expect_equal(otu_new[singleton_positions[i, 1], singleton_positions[i, 2]], 0)
    }
  }

  result_clean <- filter_occurrences_pq(data_fungi_mini, . > 1, clean_phyloseq_object = TRUE)
  
  expect_equal(sum(data_fungi_mini@otu_table)-sum(result_clean@otu_table), 
    sum(otu_orig)-sum(otu_new))
})

test_that("filter_occurrences_pq works with sample_total", {
  # Keep only values > 0.01% of sample total
  result <- filter_occurrences_pq(data_fungi, . / sample_total > 0.0001)
  expect_s4_class(result, "phyloseq")

  # Values in result should be either 0 or > 0.01% of sample total
  otu <- as(phyloseq::otu_table(result), "matrix")
  if (!phyloseq::taxa_are_rows(result)) {
    otu <- t(otu)
  }
  sample_totals <- colSums(otu)

  for (j in seq_len(ncol(otu))) {
    nonzero <- otu[, j] > 0
    if (any(nonzero) && sample_totals[j] > 0) {
      props <- otu[nonzero, j] / sample_totals[j]
      expect_true(all(props > 0.0001 | otu[nonzero, j] == 0))
    }
  }
})

test_that("mutate_occurrences_pq converts to relative abundance", {
  result <- mutate_occurrences_pq(data_fungi, . / sample_total)
  expect_s4_class(result, "phyloseq")

  otu <- as(phyloseq::otu_table(result), "matrix")
  if (!phyloseq::taxa_are_rows(result)) {
    otu <- t(otu)
  }

  # Each column should sum to 1 (or 0 if sample was empty)
  col_sums <- colSums(otu)
  expect_true(all(abs(col_sums - 1) < 1e-10 | col_sums == 0))
})

test_that("mutate_occurrences_pq applies log transformation", {
  result <- mutate_occurrences_pq(data_fungi, log1p(.))
  expect_s4_class(result, "phyloseq")

  otu_orig <- as(phyloseq::otu_table(data_fungi), "matrix")
  otu_new <- as(phyloseq::otu_table(result), "matrix")

  # Check log1p transformation was applied
  expected <- log1p(otu_orig)
  expect_equal(otu_new, expected, tolerance = 1e-10)
})

test_that("mutate_occurrences_pq handles taxa_are_rows correctly", {
  # This test ensures the function works regardless of OTU table orientation
  result <- mutate_occurrences_pq(data_fungi, . * 2)
  expect_s4_class(result, "phyloseq")

  otu_orig <- as(phyloseq::otu_table(data_fungi), "matrix")
  otu_new <- as(phyloseq::otu_table(result), "matrix")

  expect_equal(otu_new, otu_orig * 2)
  expect_equal(phyloseq::taxa_are_rows(result), phyloseq::taxa_are_rows(data_fungi))
})
