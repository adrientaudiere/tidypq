# Coerce a phyloseq tax_table into a rich data.frame

################################################################################
#' Convert a `tax_table` into a data.frame with numeric and list columns
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' The `tax_table` slot of a phyloseq object is stored as a **character**
#' matrix, so numeric information (e.g. a bootstrap support, a trait value)
#' is silently coerced to text and list-valued cells cannot be represented
#' at all. This helper lifts those limitations by returning a tibble in
#' which columns that look numeric (or logical) are converted to their
#' natural type, and selected columns can be split into genuine list
#' columns.
#'
#' It does **not** modify the phyloseq object; it produces a stand-alone
#' tibble for downstream (tidyverse) manipulation.
#'
#' @param physeq (phyloseq or taxonomyTable, required) A phyloseq object
#'   with a `tax_table` slot, or a `tax_table` object directly.
#' @param convert (logical, default TRUE) If TRUE, each column is passed
#'   through [utils::type.convert()] (with `as.is = TRUE`) so that columns
#'   holding only numbers become numeric/integer and `"TRUE"`/`"FALSE"`
#'   columns become logical. Columns with any non-numeric value are kept as
#'   character.
#' @param split (character, optional) Names of columns to turn into list
#'   columns by splitting each cell on `sep`. Useful when a rank cell packs
#'   several values (e.g. multiple GBIF matches).
#' @param sep (character, default `"/"`) Separator used to split the
#'   columns named in `split`. Passed to [base::strsplit()] with
#'   `fixed = TRUE`.
#' @param taxa_names_col (character or NULL, default `"taxon"`) Name of the
#'   leading column holding the taxa names. If NULL, taxa names are kept as
#'   row names instead of a column.
#'
#' @return A [tibble::tibble] with one row per taxon and one column per
#'   taxonomic rank, plus (unless `taxa_names_col` is NULL) a leading column
#'   of taxa names.
#' @importFrom rlang :=
#' @export
#' @author Adrien Taudière
#'
#' @seealso [phyloseq::tax_table()]
#'
#' @examples
#' data_fungi_mini <-
#'  tidypq::mutate_taxa_pq(data_fungi_mini,
#'  Mol_Abundance = taxa_sums(.))
#'
#' # Numeric-looking ranks are converted to numeric
#' df <- tax_table_to_df(data_fungi_mini)
#' df |>
#'   group_by(Class) |>
#'   summarize(mean=mean(Mol_Abundance))
#'
#'
#' \dontrun{
#' # Split a packed column into a list column
#'   tax_table_to_df(data_fungi_mini, split = "Genus", sep = "_")
#' }
tax_table_to_df <- function(
  physeq,
  convert = TRUE,
  split = NULL,
  sep = "/",
  taxa_names_col = "taxon"
) {
  if (inherits(physeq, "taxonomyTable")) {
    tt <- physeq
  } else if (methods::is(physeq, "phyloseq")) {
    tt <- phyloseq::tax_table(physeq, errorIfNULL = FALSE)
    if (is.null(tt)) {
      cli::cli_abort("{.arg physeq} has no {.field tax_table} slot.")
    }
  } else {
    cli::cli_abort(
      "{.arg physeq} must be a {.cls phyloseq} or {.cls taxonomyTable} object."
    )
  }

  taxa <- rownames(tt)
  df <- as.data.frame(tt, stringsAsFactors = FALSE)

  if (convert) {
    df[] <- lapply(df, function(col) {
      utils::type.convert(as.character(col), as.is = TRUE)
    })
  }

  if (!is.null(split)) {
    missing_cols <- setdiff(split, names(df))
    if (length(missing_cols) > 0) {
      cli::cli_abort(
        "Column{?s} {.val {missing_cols}} not found in the tax_table."
      )
    }
    for (nm in split) {
      df[[nm]] <- strsplit(as.character(df[[nm]]), sep, fixed = TRUE)
    }
  }

  out <- tibble::as_tibble(df)
  if (is.null(taxa_names_col)) {
    out <- as.data.frame(out)
    rownames(out) <- taxa
    return(out)
  }

  tibble::add_column(out, !!taxa_names_col := taxa, .before = 1)
}
