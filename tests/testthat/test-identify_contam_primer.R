skip_if_not_installed("Biostrings")

test_that("identify_contam_primer_pq returns a contam_tbl", {
  skip_if(is.null(phyloseq::refseq(data_fungi, errorIfNULL = FALSE)))
  primers <- c(fwd = "CCCTACGGGGTGCASCAG", rev = "GGACTACVSGGGTATCTAAT")
  flagged <- identify_contam_primer_pq(data_fungi, primers, verbose = FALSE)
  expect_true(is_contam_tbl(flagged))
  expect_true(all(c("matched_primers", "n_reads") %in% names(flagged)))
  if (nrow(flagged) > 0) {
    expect_true(all(flagged$method == "primer"))
  }
})

test_that("identify_contam_primer_pq flags an embedded primer", {
  skip_if(is.null(phyloseq::refseq(data_fungi, errorIfNULL = FALSE)))
  # Use a real sub-sequence of the first taxon's refseq as a fake primer
  seq1 <- as.character(phyloseq::refseq(data_fungi))[1]
  embedded <- c(planted = substr(seq1, 1, 18))
  flagged <- identify_contam_primer_pq(data_fungi, embedded, verbose = FALSE)
  expect_true(phyloseq::taxa_names(data_fungi)[1] %in% flagged$taxon)
})

test_that("identify_contam_primer_pq requires named primers and a refseq", {
  skip_if(is.null(phyloseq::refseq(data_fungi, errorIfNULL = FALSE)))
  expect_error(
    identify_contam_primer_pq(data_fungi, c("ACGT")),
    "named"
  )
  pq_no_ref <- phyloseq::phyloseq(
    phyloseq::otu_table(data_fungi),
    phyloseq::tax_table(data_fungi),
    phyloseq::sample_data(data_fungi)
  )
  expect_error(
    identify_contam_primer_pq(pq_no_ref, c(p = "ACGT")),
    "refseq"
  )
})
