skip_on_cran()

test_that("identify_contam_chimera_pq (dada2) returns a contam_tbl", {
  skip_if_not_installed("dada2")
  skip_if(is.null(phyloseq::refseq(data_fungi, errorIfNULL = FALSE)))
  flagged <- suppressMessages(
    identify_contam_chimera_pq(data_fungi, method = "dada2", verbose = FALSE)
  )
  expect_true(is_contam_tbl(flagged))
  expect_true(all(c("chimera_method", "n_reads") %in% names(flagged)))
  if (nrow(flagged) > 0) {
    expect_true(all(flagged$method == "chimera"))
    expect_true(all(flagged$chimera_method == "dada2"))
  }
})

test_that("identify_contam_chimera_pq errors without a refseq", {
  pq_no_ref <- phyloseq::phyloseq(
    phyloseq::otu_table(data_fungi),
    phyloseq::tax_table(data_fungi),
    phyloseq::sample_data(data_fungi)
  )
  expect_error(
    identify_contam_chimera_pq(pq_no_ref, method = "dada2"),
    "refseq"
  )
})

test_that("identify_contam_chimera_pq (vsearch_ref) needs an existing database", {
  skip_if(is.null(phyloseq::refseq(data_fungi, errorIfNULL = FALSE)))
  expect_error(
    identify_contam_chimera_pq(
      data_fungi,
      method = "vsearch_ref",
      database = "/nonexistent/db.fasta"
    ),
    "database"
  )
})
