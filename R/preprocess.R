#' Clean MetaPhlAn Taxonomic Names
#' 
#' @param taxa_matrix Numeric matrix of taxonomic abundances with MetaPhlAn row names.
#' @param rank Taxonomic rank to extract (default: "species").
#' @return Cleaned taxonomic matrix with standardized species names.
#' @export
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
prep_humann_pathways <- function(pathway_matrix) {
  clean_mat <- pathway_matrix[!grepl("\\|", rownames(pathway_matrix)) & 
                                !grepl("UNINTEGRATED|UNMAPPED", rownames(pathway_matrix)), , drop = FALSE]
  return(clean_mat)
}