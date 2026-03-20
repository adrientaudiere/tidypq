test_that("create_chimera_pq returns correct structure", {
  result <- create_chimera_pq(data_fungi, n_chimeras = 3)

  expect_type(result, "list")
  expect_named(result, c("physeq", "chimera_names", "parent_info", "params"))
  expect_s4_class(result$physeq, "phyloseq")
  expect_length(result$chimera_names, 3)
  expect_equal(nrow(result$parent_info), 3)
  expect_true(all(
    result$chimera_names %in% phyloseq::taxa_names(result$physeq)
  ))
})

test_that("create_chimera_pq adds chimeras to phyloseq", {
  n_orig <- phyloseq::ntaxa(data_fungi)
  result <- create_chimera_pq(data_fungi, n_chimeras = 5)

  expect_equal(phyloseq::ntaxa(result$physeq), n_orig + 5)
  expect_equal(
    phyloseq::nsamples(result$physeq),
    phyloseq::nsamples(data_fungi)
  )
})

test_that("create_chimera_pq parent_info has expected columns", {
  result <- create_chimera_pq(data_fungi, n_chimeras = 2)
  expected_cols <- c(
    "chimera",
    "parent1",
    "parent2",
    "parent_distance",
    "prop_parent1",
    "breakpoint",
    "seq_length"
  )
  expect_named(result$parent_info, expected_cols)
  expect_true(all(result$parent_info$prop_parent1 > 0))
  expect_true(all(result$parent_info$prop_parent1 < 1))
})

test_that("create_chimera_pq respects prop_mean parameter", {
  result <- create_chimera_pq(
    data_fungi,
    n_chimeras = 20,
    prop_mean = 0.7,
    prop_sd = 0.05
  )
  expect_true(mean(result$parent_info$prop_parent1) > 0.6)
})

test_that("create_chimera_pq stores params", {
  result <- create_chimera_pq(data_fungi, n_chimeras = 2, prop_mean = 0.6)
  expect_equal(result$params$prop_mean, 0.6)
  expect_equal(result$params$prop_sd, 0.15)
  expect_equal(result$params$prop_min, 0.1)
  expect_equal(result$params$min_parent_distance, 0.1)
})

test_that("create_chimera_pq works with min_parent_distance = 0", {
  result <- create_chimera_pq(
    data_fungi,
    n_chimeras = 3,
    min_parent_distance = 0
  )
  expect_s4_class(result$physeq, "phyloseq")
  expect_length(result$chimera_names, 3)
})

test_that("create_chimera_pq is reproducible with seed", {
  r1 <- create_chimera_pq(data_fungi, n_chimeras = 3, seed = 42)
  r2 <- create_chimera_pq(data_fungi, n_chimeras = 3, seed = 42)
  expect_identical(r1$chimera_names, r2$chimera_names)
  expect_identical(r1$parent_info, r2$parent_info)
})

test_that("create_chimera_pq errors without refseq", {
  pq_no_ref <- phyloseq::phyloseq(
    phyloseq::otu_table(data_fungi),
    phyloseq::tax_table(data_fungi),
    phyloseq::sample_data(data_fungi)
  )
  expect_error(create_chimera_pq(pq_no_ref), "refseq")
})

test_that("create_chimera_pq errors with invalid prop_min", {
  expect_error(
    create_chimera_pq(data_fungi, n_chimeras = 1, prop_min = 0),
    "prop_min"
  )
  expect_error(
    create_chimera_pq(data_fungi, n_chimeras = 1, prop_min = 0.5),
    "prop_min"
  )
})

test_that("create_chimera_pq errors with invalid min_parent_distance", {
  expect_error(
    create_chimera_pq(data_fungi, n_chimeras = 1, min_parent_distance = -0.1),
    "min_parent_distance"
  )
  expect_error(
    create_chimera_pq(data_fungi, n_chimeras = 1, min_parent_distance = 1),
    "min_parent_distance"
  )
})
