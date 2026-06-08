# Flag taxa whose genus belongs to a known-contaminant blocklist

################################################################################
#' Known reagent- and laboratory-contaminant genera
#'
#' @description
#' \lifecycle{experimental}
#'
#' A curated vector of bacterial genera repeatedly reported as reagent,
#' kit and laboratory contaminants in low-biomass sequencing studies. The
#' list is compiled from Salter et al. (2014), the `micRoclean` workflow
#' and the GRIMER `contaminants.yml` / `human-related.yml` resources. It is
#' intended as a starting blocklist for [contam_blocklist_pq()]; users
#' should adapt it to their own substrate and primers.
#'
#' @references
#' Salter SJ, Cox MJ, Turek EM, et al. (2014) Reagent and laboratory
#' contamination can critically impact sequence-based microbiome analyses.
#' BMC Biology 12:87. \doi{10.1186/s12915-014-0087-z}
#'
#' @return A character vector of genus names.
#'
#' @author Adrien Taudière
#' @export
#'
#' @examples
#' head(known_contaminant_genera())
#' length(known_contaminant_genera())
known_contaminant_genera <- function() {
  c(
    # Alphaproteobacteria
    "Afipia",
    "Aquabacterium",
    "Asticcacaulis",
    "Aurantimonas",
    "Beijerinckia",
    "Bosea",
    "Bradyrhizobium",
    "Brevundimonas",
    "Caulobacter",
    "Craurococcus",
    "Devosia",
    "Hoeflea",
    "Mesorhizobium",
    "Methylobacterium",
    "Novosphingobium",
    "Ochrobactrum",
    "Paracoccus",
    "Pedomicrobium",
    "Phyllobacterium",
    "Rhizobium",
    "Roseomonas",
    "Sphingobium",
    "Sphingomonas",
    "Sphingopyxis",
    # Betaproteobacteria
    "Acidovorax",
    "Azoarcus",
    "Azospira",
    "Burkholderia",
    "Comamonas",
    "Cupriavidus",
    "Curvibacter",
    "Delftia",
    "Duganella",
    "Herbaspirillum",
    "Janthinobacterium",
    "Kingella",
    "Leptothrix",
    "Limnobacter",
    "Massilia",
    "Methylophilus",
    "Methyloversatilis",
    "Oxalobacter",
    "Pelomonas",
    "Polaromonas",
    "Ralstonia",
    "Schlegelella",
    "Sulfuritalea",
    "Undibacterium",
    "Variovorax",
    # Gammaproteobacteria
    "Acinetobacter",
    "Enhydrobacter",
    "Enterobacter",
    "Escherichia",
    "Nevskia",
    "Pseudomonas",
    "Pseudoxanthomonas",
    "Psychrobacter",
    "Stenotrophomonas",
    "Xanthomonas",
    # Actinobacteria
    "Aeromicrobium",
    "Arthrobacter",
    "Brevibacterium",
    "Corynebacterium",
    "Curtobacterium",
    "Dietzia",
    "Janibacter",
    "Kocuria",
    "Microbacterium",
    "Micrococcus",
    "Propionibacterium",
    "Cutibacterium",
    "Rhodococcus",
    "Tsukamurella",
    # Firmicutes
    "Bacillus",
    "Brevibacillus",
    "Paenibacillus",
    "Staphylococcus",
    "Streptococcus",
    # Bacteroidetes
    "Chryseobacterium",
    "Dyadobacter",
    "Flavobacterium",
    "Pedobacter",
    # Deinococcus-Thermus
    "Deinococcus"
  )
}

################################################################################
#' Flag taxa whose genus belongs to a known-contaminant blocklist
#'
#' @description
#' \lifecycle{experimental}
#'
#' Flag taxa whose genus matches a blocklist of known reagent and
#' laboratory contaminants. Unlike [contam_corr_pq()] (which infers
#' contaminants from the reads-vs-depth relationship) and
#' [neg_control_classify_pq()] (which uses negative-control occurrence
#' patterns), this method relies solely on taxonomy and therefore needs no
#' negative controls. It is best used as a complementary, curation-oriented
#' filter -- review flagged taxa before removal, as several blocklist genera
#' also occur as genuine environmental taxa.
#'
#' @param physeq (phyloseq, required) A phyloseq object with a `tax_table`.
#' @param genus_blocklist (character, default:
#'   [known_contaminant_genera()]) Vector of genus names to flag.
#' @param tax_rank (character, default: "Genus") Name of the `tax_table`
#'   column holding genus-level assignments.
#' @param extra_genera (character, default: NULL) Additional genus names to
#'   append to `genus_blocklist`.
#' @param ignore_case (logical, default: TRUE) If TRUE, matching is
#'   case-insensitive.
#' @param verbose (logical, default: TRUE) If TRUE, print a `cli` summary
#'   of the flagged taxa.
#'
#' @return A named list with four elements:
#'   \describe{
#'     \item{`stats`}{A `data.frame` with one row per taxon: `taxon`,
#'       `genus`, `is_known_contaminant`, `total_reads`, `prevalence`.}
#'     \item{`contaminants`}{The subset of `stats` flagged as contaminants.}
#'     \item{`blocklist`}{The genus vector actually used.}
#'     \item{`summary`}{A named numeric vector summarizing the result.}
#'   }
#'
#' @seealso [known_contaminant_genera()], [contam_corr_pq()],
#'   [neg_control_classify_pq()]
#'
#' @author Adrien Taudière
#' @export
#'
#' @examples
#' library(MiscMetabar)
#' data(data_fungi)
#' \donttest{
#' # data_fungi is fungal, so the default bacterial blocklist matches
#' # nothing; extend it with a genus present in the data for illustration.
#' res <- contam_blocklist_pq(data_fungi, extra_genera = "Mortierella")
#' res$summary
#' head(res$contaminants)
#' }
contam_blocklist_pq <- function(
  physeq,
  genus_blocklist = known_contaminant_genera(),
  tax_rank = "Genus",
  extra_genera = NULL,
  ignore_case = TRUE,
  verbose = TRUE
) {
  MiscMetabar::verify_pq(physeq)

  if (is.null(phyloseq::tax_table(physeq, errorIfNULL = FALSE))) {
    cli::cli_abort("{.arg physeq} must contain a {.field tax_table}.")
  }
  tax <- as.data.frame(phyloseq::tax_table(physeq))
  if (!tax_rank %in% colnames(tax)) {
    cli::cli_abort(c(
      "Column {.val {tax_rank}} not found in the {.field tax_table}.",
      "i" = "Available ranks: {.val {colnames(tax)}}."
    ))
  }

  blocklist <- unique(c(genus_blocklist, extra_genera))
  blocklist <- blocklist[!is.na(blocklist) & nzchar(blocklist)]

  genus <- tax[[tax_rank]]
  if (ignore_case) {
    is_contaminant <- tolower(genus) %in% tolower(blocklist)
  } else {
    is_contaminant <- genus %in% blocklist
  }
  is_contaminant[is.na(genus)] <- FALSE

  otu <- as(phyloseq::otu_table(physeq), "matrix")
  if (!phyloseq::taxa_are_rows(physeq)) {
    otu <- t(otu)
  }
  storage.mode(otu) <- "double"
  otu <- otu[phyloseq::taxa_names(physeq), , drop = FALSE]

  stats_df <- data.frame(
    taxon = phyloseq::taxa_names(physeq),
    genus = genus,
    is_known_contaminant = is_contaminant,
    total_reads = rowSums(otu),
    prevalence = rowSums(otu > 0),
    stringsAsFactors = FALSE
  )
  stats_df <- stats_df[
    order(-stats_df$is_known_contaminant, -stats_df$total_reads),
  ]
  rownames(stats_df) <- NULL

  contaminants <- stats_df[stats_df$is_known_contaminant, , drop = FALSE]

  n_contam <- nrow(contaminants)
  n_total <- nrow(stats_df)
  summary_vec <- c(
    taxa_total = n_total,
    n_contaminants = n_contam,
    pct_contaminants = round(100 * n_contam / n_total, 1),
    reads_contaminants = sum(contaminants$total_reads),
    blocklist_size = length(blocklist)
  )

  if (verbose) {
    cli::cli_bullets(c(
      "v" = "Known-contaminant blocklist summary:",
      "*" = "Blocklist genera: {.val {length(blocklist)}}",
      "*" = "Taxa flagged: {.val {n_contam}} of {.val {n_total}} ({.val {summary_vec[['pct_contaminants']]}}%)",
      "*" = "Reads in flagged taxa: {.val {summary_vec[['reads_contaminants']]}}"
    ))
  }

  list(
    stats = stats_df,
    contaminants = contaminants,
    blocklist = blocklist,
    summary = summary_vec
  )
}
