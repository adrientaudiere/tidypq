test_that("filter_tree_pq requires phy_tree slot", {
  # Create a phyloseq without tree for testing
  physeq_no_tree <- phyloseq::phyloseq(
    phyloseq::otu_table(data_fungi),
    phyloseq::tax_table(data_fungi),
    phyloseq::sample_data(data_fungi)
  )

  expect_error(
    filter_tree_pq(physeq_no_tree),
    "does not contain a phylogenetic tree"
  )
})

# Additional tests would require a phyloseq object with a tree
# These can be added when such test data is available
