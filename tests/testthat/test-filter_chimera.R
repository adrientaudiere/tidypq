skip_on_cran()

test_that("filter_chimera_pq removes chimeric taxa in one step", {
  skip_if_not_installed("dada2")
  skip_if(is.null(phyloseq::refseq(data_fungi, errorIfNULL = FALSE)))

  flagged <- suppressMessages(
    identify_contam_chimera_pq(data_fungi, method = "dada2", verbose = FALSE)
  )
  cleaned <- suppressMessages(
    filter_chimera_pq(data_fungi, method = "dada2", verbose = FALSE)
  )

  expect_s4_class(cleaned, "phyloseq")
  expect_equal(
    phyloseq::ntaxa(cleaned),
    phyloseq::ntaxa(data_fungi) - nrow(flagged)
  )
  expect_false(any(flagged$taxon %in% phyloseq::taxa_names(cleaned)))
})

test_that("filter_chimera_pq errors without a refseq", {
  pq_no_ref <- phyloseq::phyloseq(
    phyloseq::otu_table(data_fungi),
    phyloseq::tax_table(data_fungi),
    phyloseq::sample_data(data_fungi)
  )
  expect_error(filter_chimera_pq(pq_no_ref, method = "dada2"), "refseq")
})
