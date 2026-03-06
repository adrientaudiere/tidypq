test_that("filter_taxa_pq filters by tax_table columns", {
  result <- filter_taxa_pq(data_fungi, Phylum == "Basidiomycota")
  expect_s4_class(result, "phyloseq")
  expect_true(phyloseq::ntaxa(result) < phyloseq::ntaxa(data_fungi))
  tax <- as.data.frame(phyloseq::tax_table(result))
  expect_true(all(tax$Phylum == "Basidiomycota"))
})

test_that("filter_taxa_pq filters using . pronoun", {
  min_sum <- 100
  result <- filter_taxa_pq(data_fungi, taxa_sums(.) > min_sum)
  expect_s4_class(result, "phyloseq")
  expect_true(all(phyloseq::taxa_sums(result) > min_sum))
})

test_that("filter_taxa_pq combines multiple conditions", {
  result <- filter_taxa_pq(data_fungi, Phylum == "Basidiomycota", taxa_sums(.) > 50)
  expect_s4_class(result, "phyloseq")
  tax <- as.data.frame(phyloseq::tax_table(result))
  expect_true(all(tax$Phylum == "Basidiomycota"))
  expect_true(all(phyloseq::taxa_sums(result) > 50))
})

test_that("select_taxa_pq selects columns", {
  result <- select_taxa_pq(data_fungi, Class, Phylum)
  expect_s4_class(result, "phyloseq")
  expect_equal(ncol(phyloseq::tax_table(result)), 2)
  expect_equal(colnames(phyloseq::tax_table(result)), c("Class", "Phylum"))
})

test_that("select_taxa_pq supports tidyselect range", {
  result <- select_taxa_pq(data_fungi, Phylum:Genus)
  expect_s4_class(result, "phyloseq")
  expect_true("Phylum" %in% colnames(phyloseq::tax_table(result)))
  expect_true("Class" %in% colnames(phyloseq::tax_table(result)))
  expect_true("Order" %in% colnames(phyloseq::tax_table(result)))
  expect_true("Family" %in% colnames(phyloseq::tax_table(result)))
  expect_true("Genus" %in% colnames(phyloseq::tax_table(result)))
})

test_that("mutate_taxa_pq adds new columns", {
  result <- mutate_taxa_pq(data_fungi, total_abundance = taxa_sums(.))
  expect_s4_class(result, "phyloseq")
  expect_true("total_abundance" %in% colnames(phyloseq::tax_table(result)))
  tax <- as.data.frame(phyloseq::tax_table(result))
  expected <- phyloseq::taxa_sums(data_fungi)
  expect_equal(sum(as.numeric(tax$total_abundance)- expected), 0)
})

test_that("mutate_taxa_pq modifies existing columns", {
  result <- mutate_taxa_pq(data_fungi, Genus = ifelse(is.na(Genus), "Unknown", Genus))
  expect_s4_class(result, "phyloseq")
  tax <- as.data.frame(phyloseq::tax_table(result))
  # Check that NA values were replaced
  na_count_orig <- sum(is.na(as.data.frame(phyloseq::tax_table(data_fungi))$Genus))
  na_count_new <- sum(is.na(tax$Genus))
  expect_true(na_count_new <= na_count_orig)
})

test_that("slice_taxa_pq selects by position", {
  result <- slice_taxa_pq(data_fungi, 1:10)
  expect_s4_class(result, "phyloseq")
  expect_equal(phyloseq::ntaxa(result), 10)

  original_names <- phyloseq::taxa_names(data_fungi)[1:10]
  expect_equal(phyloseq::taxa_names(result), original_names)
})

test_that("slice_taxa_pq removes by negative index", {
  n_orig <- phyloseq::ntaxa(data_fungi)
  result <- slice_taxa_pq(data_fungi, -(1:5))
  expect_equal(phyloseq::ntaxa(result), n_orig - 5)
})

test_that("arrange_taxa_pq reorders taxa by taxonomy", {
  result <- arrange_taxa_pq(data_fungi, Phylum)
  expect_s4_class(result, "phyloseq")
  tax <- as.data.frame(phyloseq::tax_table(result))
  expect_true(!is.unsorted(tax$Phylum, na.rm = TRUE))
})

test_that("arrange_taxa_pq works with . pronoun", {
  result <- arrange_taxa_pq(data_fungi, taxa_sums(.))
  expect_s4_class(result, "phyloseq")
  sums <- phyloseq::taxa_sums(result)
  expect_true(!is.unsorted(sums))
})

test_that("rename_taxa_pq renames columns", {
  result <- rename_taxa_pq(data_fungi, tax_Phylum = Phylum)
  expect_s4_class(result, "phyloseq")
  expect_true("tax_Phylum" %in% colnames(phyloseq::tax_table(result)))
  expect_false("Phylum" %in% colnames(phyloseq::tax_table(result)))
})

test_that("arrange_taxa_pq works with phy_tree slot", {
  # Use a subset of GlobalPatterns which has a phy_tree
  data("GlobalPatterns", package = "phyloseq")
  gp_small <- phyloseq::prune_taxa(phyloseq::taxa_names(GlobalPatterns)[1:50], GlobalPatterns)

  # Verify the test object has a phy_tree

  expect_false(is.null(phyloseq::phy_tree(gp_small, errorIfNULL = FALSE)))

  result <- arrange_taxa_pq(gp_small, Genus)

  expect_s4_class(result, "phyloseq")
  expect_false(is.null(phyloseq::phy_tree(result, errorIfNULL = FALSE)))

  # Check taxa order matches in tax_table and otu_table
  expect_equal(phyloseq::taxa_names(result), rownames(phyloseq::tax_table(result)))

  # Check phy_tree tips match taxa names (order may differ due to tree topology)
  expect_setequal(phyloseq::phy_tree(result)$tip.label, phyloseq::taxa_names(result))

  # Check that tax_table is sorted by Phylum
  tax <- as.data.frame(phyloseq::tax_table(result))
  expect_true(!is.unsorted(tax$Phylum, na.rm = TRUE))
})
