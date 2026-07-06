test_that("tax_table_to_df returns a tibble with a taxa column", {
  df <- tax_table_to_df(data_fungi_mini)
  expect_s3_class(df, "tbl_df")
  expect_equal(nrow(df), phyloseq::ntaxa(data_fungi_mini))
  expect_equal(names(df)[1], "taxon")
  expect_setequal(df$taxon, phyloseq::taxa_names(data_fungi_mini))
})

test_that("tax_table_to_df keeps taxa as row names when taxa_names_col is NULL", {
  df <- tax_table_to_df(data_fungi_mini, taxa_names_col = NULL)
  expect_s3_class(df, "data.frame")
  expect_setequal(rownames(df), phyloseq::taxa_names(data_fungi_mini))
  expect_false("taxon" %in% names(df))
})

test_that("tax_table_to_df converts numeric-looking columns to numeric", {
  pq <- data_fungi_mini
  tt <- as.data.frame(phyloseq::tax_table(pq))
  tt$Support <- as.character(seq_len(nrow(tt)))
  phyloseq::tax_table(pq) <- phyloseq::tax_table(as.matrix(tt))

  converted <- tax_table_to_df(pq)
  expect_type(converted$Support, "integer")

  not_converted <- tax_table_to_df(pq, convert = FALSE)
  expect_type(not_converted$Support, "character")
})

test_that("tax_table_to_df can split a column into a list column", {
  pq <- data_fungi_mini
  tt <- as.data.frame(phyloseq::tax_table(pq))
  tt$Multi <- paste(tt[[1]], tt[[2]], sep = "/")
  phyloseq::tax_table(pq) <- phyloseq::tax_table(as.matrix(tt))

  df <- tax_table_to_df(pq, split = "Multi")
  expect_type(df$Multi, "list")
  expect_true(all(lengths(df$Multi) == 2))
})

test_that("tax_table_to_df accepts a tax_table object directly", {
  tt <- phyloseq::tax_table(data_fungi_mini)
  df <- tax_table_to_df(tt)
  expect_s3_class(df, "tbl_df")
  expect_equal(nrow(df), phyloseq::ntaxa(data_fungi_mini))
})

test_that("tax_table_to_df errors on bad input and missing split column", {
  expect_error(tax_table_to_df(1:10), "phyloseq")
  expect_error(
    tax_table_to_df(data_fungi_mini, split = "not_a_column"),
    "not found"
  )
})
