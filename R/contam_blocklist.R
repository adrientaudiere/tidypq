# Curated blocklist of known reagent and laboratory contaminant genera

################################################################################
#' Known reagent- and laboratory-contaminant genera
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' A curated vector of bacterial genera repeatedly reported as reagent,
#' kit and laboratory contaminants in low-biomass sequencing studies. The
#' list is compiled from Salter et al. (2014), the `micRoclean` workflow
#' and the GRIMER `contaminants.yml` / `human-related.yml` resources. It is
#' intended as a starting blocklist for [identify_contam_blocklist_pq()];
#' users should adapt it to their own substrate and primers.
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
