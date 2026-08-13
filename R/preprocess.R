#' Clean MetaPhlAn Taxonomic Names
#'
#' @param taxa_matrix Numeric matrix of taxonomic abundances with MetaPhlAn row names.
#' @param rank Taxonomic rank to extract (default: "species").
#' @return Cleaned taxonomic matrix with standardized species names.
#' @export
#'
#' @examples
#' dummy_taxa <- matrix(
#'   c(10.5, 89.5, 20.0, 80.0),
#'   nrow = 2, ncol = 2,
#'   dimnames = list(
#'     c(
#'       "k__Bacteria|s__Staphylococcus_aureus",
#'       "k__Bacteria|s__Bacteroides_fragilis"
#'     ),
#'     c("Sample1", "Sample2")
#'   )
#' )
#' prep_metaphlan_taxa(dummy_taxa, rank = "species")
prep_metaphlan_taxa <- function(taxa_matrix, rank = "species") {
  prefix <- paste0(substr(rank, 1, 1), "__")
  clean_mat <- taxa_matrix[grepl(prefix, rownames(taxa_matrix)) & !grepl("t__", rownames(taxa_matrix)), , drop = FALSE]
  rownames(clean_mat) <- chartr("_", " ", sub(paste0(".*", prefix), "", rownames(clean_mat)))
  return(clean_mat)
}

#' Clean HUMAnN Pathway Matrix
#'
#' @param pathway_matrix Numeric matrix of pathway abundances from HUMAnN.
#' @return Cleaned pathway matrix containing only community-level unstratified pathways.
#' @export
#'
#' @examples
#' dummy_pathways <- matrix(
#'   c(150.2, 300.8, 0, 120.0, 450.5, 10.0),
#'   nrow = 3, ncol = 2,
#'   dimnames = list(
#'     c("PWY-101: gnb1", "UNMAPPED", "PWY-101|s__Staphylococcus_aureus"),
#'     c("Sample1", "Sample2")
#'   )
#' )
#' prep_humann_pathways(dummy_pathways)
prep_humann_pathways <- function(pathway_matrix) {
  clean_mat <- pathway_matrix[!grepl("\\|", rownames(pathway_matrix)) &
    !grepl("UNINTEGRATED|UNMAPPED", rownames(pathway_matrix)), , drop = FALSE]
  return(clean_mat)
}
