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
  result <- filter_taxa_pq(
    data_fungi,
    Phylum == "Basidiomycota",
    taxa_sums(.) > 50
  )
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

test_that("select_taxa_pq errors with no arguments", {
  expect_error(
    select_taxa_pq(data_fungi),
    "No columns selected"
  )
})

test_that("mutate_taxa_pq adds new columns", {
  result <- mutate_taxa_pq(data_fungi, total_abundance = taxa_sums(.))
  expect_s4_class(result, "phyloseq")
  expect_true("total_abundance" %in% colnames(phyloseq::tax_table(result)))
  tax <- as.data.frame(phyloseq::tax_table(result))
  expected <- phyloseq::taxa_sums(data_fungi)
  expect_equal(as.numeric(tax$total_abundance), as.numeric(expected))
})

test_that("mutate_taxa_pq modifies existing columns", {
  result <- mutate_taxa_pq(
    data_fungi,
    Genus = ifelse(is.na(Genus), "Unknown", Genus)
  )
  expect_s4_class(result, "phyloseq")
  tax <- as.data.frame(phyloseq::tax_table(result))
  # Check that NA values were replaced
  na_count_orig <- sum(is.na(
    as.data.frame(phyloseq::tax_table(data_fungi))$Genus
  ))
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

test_that("arrange_taxa_pq preserves tree topology after reordering", {
  # Use a subset of GlobalPatterns which has a phy_tree
  data("GlobalPatterns", package = "phyloseq")
  gp_small <- phyloseq::prune_taxa(
    phyloseq::taxa_names(GlobalPatterns)[1:50],
    GlobalPatterns
  )

  dist_before <- ape::cophenetic.phylo(phyloseq::phy_tree(gp_small))
  result <- arrange_taxa_pq(gp_small, Genus)
  dist_after <- ape::cophenetic.phylo(phyloseq::phy_tree(result))

  # Same taxa
  expect_setequal(rownames(dist_after), rownames(dist_before))

  # Pairwise distances preserved (reorder both matrices to same taxa order)
  taxa_order <- rownames(dist_before)
  expect_equal(
    dist_after[taxa_order, taxa_order],
    dist_before[taxa_order, taxa_order]
  )
})

test_that("arrange_taxa_pq works with phy_tree slot", {
  # Use a subset of GlobalPatterns which has a phy_tree
  data("GlobalPatterns", package = "phyloseq")
  gp_small <- phyloseq::prune_taxa(
    phyloseq::taxa_names(GlobalPatterns)[1:50],
    GlobalPatterns
  )

  # Verify the test object has a phy_tree

  expect_false(is.null(phyloseq::phy_tree(gp_small, errorIfNULL = FALSE)))

  result <- arrange_taxa_pq(gp_small, Genus)

  expect_s4_class(result, "phyloseq")
  expect_false(is.null(phyloseq::phy_tree(result, errorIfNULL = FALSE)))

  # Check taxa order matches in tax_table and otu_table
  expect_equal(
    phyloseq::taxa_names(result),
    rownames(phyloseq::tax_table(result))
  )

  # Check phy_tree tips match taxa names (order may differ due to tree topology)
  expect_setequal(
    phyloseq::phy_tree(result)$tip.label,
    phyloseq::taxa_names(result)
  )

  # Check that tax_table is sorted by Phylum
  tax <- as.data.frame(phyloseq::tax_table(result))
  expect_true(!is.unsorted(tax$Phylum, na.rm = TRUE))
})

# ---- Regression: phyloseq objects with desynced otu_table/tax_table order ----
# Such objects are produced e.g. by MiscMetabar::mumu_pq(), which assigns the
# otu_table directly in mumu output order while the tax_table stays in input
# order. taxa_names() prioritises the otu_table, so verbs must not assume that
# tax_table rows are positionally aligned with taxa_names().

make_desynced_pq <- function() {
  taxa <- paste0("ASV", 1:6)
  tax <- data.frame(
    Kingdom = c("Archaea", "Bacteria", "Archaea", "Bacteria", "Archaea", "Bacteria"),
    Phylum = c("Thaum", "Firm", "Eury", "Bact", "Cren", "Bact"),
    row.names = taxa,
    stringsAsFactors = FALSE
  )
  otu <- matrix(
    c(10, 0, 5, 0, 3, 0,
      0, 7, 0, 2, 0, 4,
      4, 0, 8, 0, 1, 0),
    nrow = 6, byrow = FALSE,
    dimnames = list(taxa, c("S1", "S2", "S3"))
  )
  pq <- phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = TRUE),
    phyloseq::tax_table(as.matrix(tax))
  )
  # Desync: tax_table in a different order than the otu_table (same taxa set).
  shuffled <- c("ASV3", "ASV1", "ASV5", "ASV2", "ASV4", "ASV6")
  pq@tax_table <- phyloseq::tax_table(as.matrix(tax[shuffled, , drop = FALSE]))
  pq
}

test_that("filter_taxa_pq is correct on a desynced object", {
  pq <- make_desynced_pq()
  expect_false(identical(
    phyloseq::taxa_names(pq),
    rownames(phyloseq::tax_table(pq))
  ))
  result <- filter_taxa_pq(pq, Kingdom == "Archaea", clean_phyloseq_object = FALSE)
  expect_s4_class(result, "phyloseq")
  kingdoms <- as.character(phyloseq::tax_table(result)[, "Kingdom"])
  expect_equal(unique(kingdoms), "Archaea")
  expect_setequal(phyloseq::taxa_names(result), c("ASV1", "ASV3", "ASV5"))
})

test_that("mutate_taxa_pq is correct on a desynced object", {
  pq <- make_desynced_pq()
  result <- mutate_taxa_pq(pq, total = phyloseq::taxa_sums(.))
  tax <- as.data.frame(phyloseq::tax_table(result))
  expected <- phyloseq::taxa_sums(pq)
  expect_equal(as.numeric(tax$total), as.numeric(expected))
  # Values must be aligned to the correct taxon, not the tax_table order
  expect_equal(
    as.numeric(tax[phyloseq::taxa_names(pq), "total"]),
    as.numeric(expected[phyloseq::taxa_names(pq)])
  )
})

test_that("slice_taxa_pq is correct on a desynced object", {
  pq <- make_desynced_pq()
  result <- slice_taxa_pq(pq, 1:3, clean_phyloseq_object = FALSE)
  # Positions refer to taxa_names() order, not the (desynced) tax_table order
  expect_equal(
    phyloseq::taxa_names(result),
    phyloseq::taxa_names(pq)[1:3]
  )
})

test_that("arrange_taxa_pq is correct on a desynced object", {
  pq <- make_desynced_pq()
  result <- arrange_taxa_pq(pq, Kingdom, clean_phyloseq_object = FALSE)
  tax <- as.data.frame(phyloseq::tax_table(result))
  expect_true(!is.unsorted(tax$Kingdom, na.rm = TRUE))
  # tax_table and otu_table must stay aligned after arranging
  expect_equal(
    phyloseq::taxa_names(result),
    rownames(phyloseq::tax_table(result))
  )
})
