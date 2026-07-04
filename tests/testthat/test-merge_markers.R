make_marker_pq <- function(taxa, ranks, samples, seed = 1) {
  set.seed(seed)
  otu <- matrix(
    sample(0:20, length(taxa) * length(samples), replace = TRUE),
    nrow = length(taxa),
    dimnames = list(taxa, samples)
  )
  tax <- matrix(
    paste(rep(ranks, each = length(taxa)), seq_along(taxa), sep = "_"),
    nrow = length(taxa),
    dimnames = list(taxa, ranks)
  )
  sam <- data.frame(
    row.names = samples,
    site = rep_len(c("A", "B"), length(samples)),
    depth = seq_along(samples)
  )
  phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = TRUE),
    phyloseq::tax_table(tax),
    phyloseq::sample_data(sam)
  )
}

test_that("merge_markers_pq disambiguates taxa and never merges them", {
  p1 <- make_marker_pq(
    c("ASV_1", "ASV_2"),
    c("Kingdom", "Phylum"),
    c("s1", "s2", "s3")
  )
  p2 <- make_marker_pq(
    c("ASV_1", "ASV_3"),
    c("Kingdom", "Phylum", "Genus"),
    c("s1", "s2", "s3"),
    seed = 2
  )

  m <- suppressMessages(merge_markers_pq(list(m16S = p1, mITS = p2)))

  expect_equal(phyloseq::ntaxa(m), 4)
  expect_setequal(
    phyloseq::taxa_names(m),
    c("m16S_ASV_1", "m16S_ASV_2", "mITS_ASV_1", "mITS_ASV_3")
  )
  expect_equal(phyloseq::nsamples(m), 3)
})

test_that("merge_markers_pq records the marker of origin in tax_table", {
  p1 <- make_marker_pq(c("ASV_1"), c("Kingdom"), c("s1", "s2"))
  p2 <- make_marker_pq(c("ASV_1"), c("Kingdom"), c("s1", "s2"), seed = 2)

  m <- suppressMessages(merge_markers_pq(list(m16S = p1, mITS = p2)))
  markers <- as.vector(phyloseq::tax_table(m)[, "marker"])

  expect_setequal(markers, c("m16S", "mITS"))
  expect_true("marker" %in% colnames(phyloseq::tax_table(m)))
})

test_that("merge_markers_pq unions samples and zero-fills unobserved occurrences", {
  p1 <- make_marker_pq(c("ASV_1"), c("Kingdom"), c("s1", "s2"))
  p2 <- make_marker_pq(c("ASV_9"), c("Kingdom"), c("s2", "s3"), seed = 2)

  m <- suppressMessages(merge_markers_pq(list(a = p1, b = p2)))
  otu <- as(phyloseq::otu_table(m), "matrix")

  expect_setequal(colnames(otu), c("s1", "s2", "s3"))
  # a's taxon is not observable in s3 (only sequenced by b) -> 0
  expect_equal(otu["a_ASV_1", "s3"], 0)
  # b's taxon is not observable in s1 (only sequenced by a) -> 0
  expect_equal(otu["b_ASV_9", "s1"], 0)
})

test_that("merge_markers_pq full-joins sam_data columns", {
  p1 <- make_marker_pq(c("ASV_1"), c("Kingdom"), c("s1", "s2"))
  p2 <- make_marker_pq(c("ASV_9"), c("Kingdom"), c("s2", "s3"), seed = 2)
  phyloseq::sample_data(p2)$extra <- c("u", "v")

  m <- suppressMessages(merge_markers_pq(list(a = p1, b = p2)))
  sam <- as(phyloseq::sample_data(m), "data.frame")

  expect_true(all(c("site", "depth", "extra") %in% colnames(sam)))
  expect_equal(nrow(sam), 3)
  # s1 has no `extra` value (absent from marker b) -> NA
  expect_true(is.na(sam["s1", "extra"]))
})

test_that("merge_markers_pq accepts marker_names for an unnamed list", {
  p1 <- make_marker_pq(c("ASV_1"), c("Kingdom"), c("s1", "s2"))
  p2 <- make_marker_pq(c("ASV_1"), c("Kingdom"), c("s1", "s2"), seed = 2)

  m <- suppressMessages(
    merge_markers_pq(list(p1, p2), marker_names = c("16S", "ITS"))
  )
  expect_setequal(
    phyloseq::taxa_names(m),
    c("16S_ASV_1", "ITS_ASV_1")
  )
})

test_that("merge_markers_pq errors on colliding names when disambiguate is off", {
  p1 <- make_marker_pq(c("ASV_1"), c("Kingdom"), c("s1", "s2"))
  p2 <- make_marker_pq(c("ASV_1"), c("Kingdom"), c("s1", "s2"), seed = 2)

  expect_error(
    suppressMessages(
      merge_markers_pq(list(a = p1, b = p2), disambiguate_taxa = FALSE)
    ),
    "collide"
  )
})

test_that("merge_markers_pq validates its inputs", {
  p1 <- make_marker_pq(c("ASV_1"), c("Kingdom"), c("s1", "s2"))

  expect_error(merge_markers_pq(list(p1)), "at least two")
  expect_error(merge_markers_pq(list(p1, "not_pq")), "phyloseq objects")
  expect_error(
    merge_markers_pq(list(p1, p1)),
    "marker name"
  )
  expect_error(
    merge_markers_pq(list(a = p1, b = p1), marker_names = c("x", "x")),
    "unique"
  )
})
