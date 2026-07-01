data(data_fungi, package = "MiscMetabar")

make_small <- function() {
  pq <- phyloseq::prune_samples(
    phyloseq::sample_names(data_fungi)[1:6],
    data_fungi
  )
  phyloseq::prune_taxa(phyloseq::taxa_names(pq)[1:8], pq)
}

test_that("canonicalize_pq_order reorders desynchronised sam_data with a message", {
  pq <- make_small()
  pq@sam_data <- pq@sam_data[rev(phyloseq::sample_names(pq)), ]

  expect_message(
    res <- canonicalize_pq_order(pq),
    "Reordering sam_data rows"
  )
  expect_identical(
    phyloseq::sample_names(res@otu_table),
    rownames(res@sam_data)
  )
})

test_that("canonicalize_pq_order reorders desynchronised tax_table with a message", {
  pq <- make_small()
  pq@tax_table <- pq@tax_table[rev(phyloseq::taxa_names(pq)), ]

  expect_message(
    res <- canonicalize_pq_order(pq),
    "Reordering tax_table rows"
  )
  expect_identical(
    phyloseq::taxa_names(res@otu_table),
    rownames(res@tax_table)
  )
})

test_that("canonicalize_pq_order is a silent no-op on a consistent object", {
  pq <- make_small()
  expect_no_message(res <- canonicalize_pq_order(pq))
  expect_identical(res, pq)
})

test_that("mutators self-heal a misordered object and emit a consistent result", {
  pq <- make_small()
  pq@sam_data <- pq@sam_data[rev(phyloseq::sample_names(pq)), ]

  res <- suppressMessages(mutate_samdata_pq(
    pq,
    foo = seq_len(phyloseq::nsamples(pq))
  ))
  expect_silent(MiscMetabar::verify_pq(res))

  pq2 <- make_small()
  pq2@tax_table <- pq2@tax_table[rev(phyloseq::taxa_names(pq2)), ]
  res2 <- suppressMessages(mutate_taxa_pq(pq2, newrank = "x"))
  expect_silent(MiscMetabar::verify_pq(res2))
})
