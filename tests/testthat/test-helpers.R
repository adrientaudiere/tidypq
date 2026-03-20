test_that("taxa_prevalence returns correct prevalence", {
  prev <- taxa_prevalence(data_fungi)
  expect_type(prev, "integer")
  expect_length(prev, phyloseq::ntaxa(data_fungi))
  expect_true(all(names(prev) == phyloseq::taxa_names(data_fungi)))

  # Prevalence should be between 0 and number of samples
  expect_true(all(prev >= 0))
  expect_true(all(prev <= phyloseq::nsamples(data_fungi)))
})

test_that("taxa_prevalence respects threshold", {
  prev0 <- taxa_prevalence(data_fungi, threshold = 0)
  prev10 <- taxa_prevalence(data_fungi, threshold = 10)

  # Higher threshold should give equal or lower prevalence
  expect_true(all(prev10 <= prev0))
})

test_that("taxa_prevalence works with filter_taxa_pq", {
  # Keep taxa present in at least 5 samples
  result <- filter_taxa_pq(data_fungi, taxa_prevalence(.) >= 5)
  expect_s4_class(result, "phyloseq")

  # Verify all remaining taxa have prevalence >= 5
  prev <- taxa_prevalence(result)
  expect_true(all(prev >= 5))
})

################################################################################
# decontam_sam_control tests
################################################################################

test_that("decontam_sam_control removes contamination based on controls", {
  # Add control column: first 2 samples are controls
  pq <- mutate_samdata_pq(data_fungi, is_control = seq_len(nsamples(.)) <= 2)

  result <- decontam_sam_control(pq, is_control, verbose = FALSE)
  expect_s4_class(result, "phyloseq")

  # Control samples should be kept by default (remove_controls = FALSE)
  expect_equal(phyloseq::nsamples(result), phyloseq::nsamples(data_fungi))
})

test_that("decontam_sam_control removes controls when remove_controls = TRUE", {
  pq <- mutate_samdata_pq(data_fungi, is_control = seq_len(nsamples(.)) <= 2)

  result <- decontam_sam_control(
    pq,
    is_control,
    remove_controls = TRUE,
    verbose = FALSE
  )
  expect_s4_class(result, "phyloseq")

  # Control samples should be removed
  expect_equal(phyloseq::nsamples(result), phyloseq::nsamples(data_fungi) - 2)
  expect_false(any(
    phyloseq::sample_names(result) %in% phyloseq::sample_names(pq)[1:2]
  ))
})

test_that("decontam_sam_control works with different summary functions", {
  pq <- mutate_samdata_pq(data_fungi, is_control = seq_len(nsamples(.)) <= 2)

  # max (default) - most conservative
  result_max <- decontam_sam_control(pq, is_control, fun = max, verbose = FALSE)

  # mean - less conservative
  result_mean <- decontam_sam_control(
    pq,
    is_control,
    fun = mean,
    verbose = FALSE
  )

  expect_s4_class(result_max, "phyloseq")
  expect_s4_class(result_mean, "phyloseq")

  # max threshold should remove at least as many occurrences as mean
  # (max >= mean for positive values)
  expect_true(
    sum(phyloseq::otu_table(result_max)) <=
      sum(phyloseq::otu_table(result_mean))
  )
})

test_that("decontam_sam_control works with global_threshold = TRUE", {
  pq <- mutate_samdata_pq(data_fungi, is_control = seq_len(nsamples(.)) <= 2)

  # Per-taxon thresholds (default)
  result_per_taxon <- decontam_sam_control(
    pq,
    is_control,
    global_threshold = FALSE,
    verbose = FALSE
  )

  # Global threshold
  result_global <- decontam_sam_control(
    pq,
    is_control,
    global_threshold = TRUE,
    verbose = FALSE
  )

  expect_s4_class(result_per_taxon, "phyloseq")
  expect_s4_class(result_global, "phyloseq")

  # Both should produce valid phyloseq objects
  expect_true(phyloseq::nsamples(result_global) > 0)
})

test_that("decontam_sam_control works with custom function", {
  pq <- mutate_samdata_pq(data_fungi, is_control = seq_len(nsamples(.)) <= 2)

  # Custom function: 2x the max
  result <- decontam_sam_control(
    pq,
    is_control,
    fun = \(x) 2 * max(x),
    verbose = FALSE
  )
  expect_s4_class(result, "phyloseq")
})

test_that("decontam_sam_control errors with no control samples", {
  pq <- mutate_samdata_pq(data_fungi, is_control = FALSE)

  expect_error(
    decontam_sam_control(pq, is_control, verbose = FALSE),
    "No samples match"
  )
})

test_that("decontam_sam_control errors with all control samples", {
  pq <- mutate_samdata_pq(data_fungi, is_control = TRUE)

  expect_error(
    decontam_sam_control(pq, is_control, verbose = FALSE),
    "All samples are controls"
  )
})

test_that("decontam_sam_control errors with non-logical condition", {
  pq <- mutate_samdata_pq(data_fungi, control_value = 1)

  expect_error(
    decontam_sam_control(pq, control_value, verbose = FALSE),
    "must evaluate to a logical"
  )
})

test_that("decontam_sam_control correctly sets values to zero", {
  # Create a simple controlled test case
  pq <- mutate_samdata_pq(data_fungi, is_control = seq_len(nsamples(.)) <= 2)

  # Get the control threshold for the first taxon (using max)
  tar <- phyloseq::taxa_are_rows(pq)
  otu <- as(phyloseq::otu_table(pq), "matrix")
  if (!tar) {
    otu <- t(otu)
  }

  control_samples <- phyloseq::sample_names(pq)[1:2]
  first_taxon <- phyloseq::taxa_names(pq)[1]
  threshold <- max(otu[first_taxon, control_samples])

  result <- decontam_sam_control(
    pq,
    is_control,
    clean_phyloseq_object = FALSE,
    verbose = FALSE
  )

  # Get result OTU table
  result_otu <- as(phyloseq::otu_table(result), "matrix")
  if (!phyloseq::taxa_are_rows(result)) {
    result_otu <- t(result_otu)
  }

  # For non-control samples, values <= threshold should now be 0
  # (If the taxon still exists in result - it might have been removed)
  if (first_taxon %in% rownames(result_otu)) {
    non_control <- setdiff(colnames(result_otu), control_samples)
    remaining_values <- result_otu[first_taxon, non_control]
    expect_true(all(remaining_values == 0 | remaining_values > threshold))
  }
})

################################################################################
# decontam_taxa_control tests
################################################################################

test_that("decontam_taxa_control removes contamination based on control taxa", {
  # Use first 3 taxa as controls via taxa_names
  control_taxa <- phyloseq::taxa_names(data_fungi)[1:3]

  result <- decontam_taxa_control(
    data_fungi,
    taxa_names(.) %in% control_taxa,
    verbose = FALSE,
    clean_phyloseq_object = FALSE
  )
  expect_s4_class(result, "phyloseq")

  # Control taxa should be removed by default
  expect_equal(phyloseq::ntaxa(result), phyloseq::ntaxa(data_fungi) - 3)
  expect_false(any(control_taxa %in% phyloseq::taxa_names(result)))
})

test_that("decontam_taxa_control works with tax_table condition", {
  # Use a condition on tax_table (e.g., select by Genus)
  # First find a Genus that exists
  tax_mat <- as(phyloseq::tax_table(data_fungi), "matrix")
  test_genus <- tax_mat[1, "Genus"]

  result <- decontam_taxa_control(
    data_fungi,
    Genus == test_genus,
    verbose = FALSE
  )
  expect_s4_class(result, "phyloseq")
})

test_that("decontam_taxa_control keeps control taxa when remove_control_taxa = FALSE", {
  control_taxa <- phyloseq::taxa_names(data_fungi)[1:3]

  result <- decontam_taxa_control(
    data_fungi,
    taxa_names(.) %in% control_taxa,
    remove_control_taxa = FALSE,
    verbose = FALSE
  )
  expect_s4_class(result, "phyloseq")

  # Control taxa should still be present
  expect_true(all(control_taxa %in% phyloseq::taxa_names(result)))
})

test_that("decontam_taxa_control works with global_threshold = TRUE", {
  control_taxa <- phyloseq::taxa_names(data_fungi)[1:3]

  # Per-sample thresholds (default)
  result_per_sample <- decontam_taxa_control(
    data_fungi,
    taxa_names(.) %in% control_taxa,
    global_threshold = FALSE,
    verbose = FALSE
  )

  # Global threshold
  result_global <- decontam_taxa_control(
    data_fungi,
    taxa_names(.) %in% control_taxa,
    global_threshold = TRUE,
    verbose = FALSE
  )

  expect_s4_class(result_per_sample, "phyloseq")
  expect_s4_class(result_global, "phyloseq")
})

test_that("decontam_taxa_control works with different summary functions", {
  control_taxa <- phyloseq::taxa_names(data_fungi)[1:3]

  result_max <- decontam_taxa_control(
    data_fungi,
    taxa_names(.) %in% control_taxa,
    fun = max,
    verbose = FALSE
  )
  result_mean <- decontam_taxa_control(
    data_fungi,
    taxa_names(.) %in% control_taxa,
    fun = mean,
    verbose = FALSE
  )

  expect_s4_class(result_max, "phyloseq")
  expect_s4_class(result_mean, "phyloseq")

  # max threshold should remove at least as many occurrences as mean
  expect_true(
    sum(phyloseq::otu_table(result_max)) <=
      sum(phyloseq::otu_table(result_mean))
  )
})

test_that("decontam_taxa_control errors with no matching taxa", {
  expect_error(
    decontam_taxa_control(
      data_fungi,
      Genus == "NONEXISTENT_GENUS_XYZ",
      verbose = FALSE
    ),
    "No taxa match"
  )
})

test_that("decontam_taxa_control errors with all control taxa", {
  expect_error(
    decontam_taxa_control(
      data_fungi,
      taxa_names(.) %in% taxa_names(.),
      verbose = FALSE
    ),
    "All taxa are control taxa"
  )
})

test_that("decontam_taxa_control errors with non-logical condition", {
  expect_error(
    decontam_taxa_control(data_fungi, taxa_sums(.), verbose = FALSE),
    "must evaluate to a logical"
  )
})

test_that("decontam_sam_control verbose output works", {
  pq <- mutate_samdata_pq(
    data_fungi,
    is_control = seq_len(nsamples(.)) <= 2
  )
  expect_message(
    decontam_sam_control(pq, is_control, verbose = TRUE),
    "Decontamination complete"
  )
})

test_that("decontam_sam_control verbose with global_threshold", {
  pq <- mutate_samdata_pq(
    data_fungi,
    is_control = seq_len(nsamples(.)) <= 2
  )
  expect_message(
    decontam_sam_control(
      pq,
      is_control,
      global_threshold = TRUE,
      verbose = TRUE
    ),
    "Global threshold value"
  )
})

test_that("decontam_taxa_control verbose output works", {
  control_taxa <- phyloseq::taxa_names(data_fungi)[1:3]
  expect_message(
    decontam_taxa_control(
      data_fungi,
      taxa_names(.) %in% control_taxa,
      verbose = TRUE
    ),
    "Decontamination complete"
  )
})

test_that("decontam_taxa_control verbose with global_threshold", {
  control_taxa <- phyloseq::taxa_names(data_fungi)[1:3]
  expect_message(
    decontam_taxa_control(
      data_fungi,
      taxa_names(.) %in% control_taxa,
      global_threshold = TRUE,
      verbose = TRUE
    ),
    "Global threshold value"
  )
})
