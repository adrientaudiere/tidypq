test_that("filter_tree_pq requires phy_tree slot", {
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

# Helper to build a phyloseq with tree for these tests
make_pq_with_tree <- function() {
  taxa_nms <- phyloseq::taxa_names(data_fungi)
  tree <- ape::rtree(length(taxa_nms), tip.label = taxa_nms)
  phyloseq::phyloseq(
    phyloseq::otu_table(data_fungi),
    phyloseq::tax_table(data_fungi),
    phyloseq::sample_data(data_fungi),
    phyloseq::phy_tree(tree)
  )
}

test_that("filter_tree_pq synchronizes tree with OTU table", {
  pq <- make_pq_with_tree()
  result <- filter_tree_pq(pq)

  expect_s4_class(result, "phyloseq")
  tips <- result@phy_tree$tip.label
  expect_true(all(tips %in% phyloseq::taxa_names(result)))
})

test_that("filter_tree_pq filters by explicit taxa", {
  pq <- make_pq_with_tree()
  keep <- phyloseq::taxa_names(pq)[1:5]

  result <- filter_tree_pq(pq, taxa = keep)

  expect_equal(phyloseq::ntaxa(result), 5)
  expect_true(all(phyloseq::taxa_names(result) %in% keep))
  expect_equal(length(result@phy_tree$tip.label), 5)
})

test_that("filter_tree_pq filters by pattern", {
  pq <- make_pq_with_tree()
  first_taxon <- phyloseq::taxa_names(pq)[1]
  pat <- substr(first_taxon, 1, 5)

  result <- filter_tree_pq(pq, pattern = pat)
  expect_s4_class(result, "phyloseq")
  expect_true(all(grepl(pat, phyloseq::taxa_names(result))))
})

test_that("filter_tree_pq pattern with invert = TRUE", {
  pq <- make_pq_with_tree()
  first_taxon <- phyloseq::taxa_names(pq)[1]
  pat <- paste0("^", first_taxon, "$")

  result <- filter_tree_pq(pq, pattern = pat, invert = TRUE)
  expect_false(first_taxon %in% phyloseq::taxa_names(result))
})

test_that("filter_tree_pq errors when no taxa remain", {
  pq <- make_pq_with_tree()
  expect_error(
    filter_tree_pq(pq, taxa = "NONEXISTENT"),
    "No taxa remain"
  )
})

test_that("filter_tree_pq respects clean_phyloseq_object", {
  pq <- make_pq_with_tree()
  keep <- phyloseq::taxa_names(pq)[1:3]

  result <- filter_tree_pq(
    pq,
    taxa = keep,
    clean_phyloseq_object = FALSE
  )
  expect_s4_class(result, "phyloseq")
  expect_equal(length(result@phy_tree$tip.label), 3)
})
