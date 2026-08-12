#' Synthetic Example Taxa Abundance Matrix
#'
#' A simulated count matrix representing bacterial taxonomic abundances across samples.
#'
#' @format A numeric matrix with 20 rows (taxa) and 30 columns (samples).
#' @source Simulated synthetic data for MiCARA testing and demonstration.
"example_taxa"

#' Synthetic Example Functional Pathways Matrix
#'
#' A simulated count matrix representing functional pathway abundances across samples.
#'
#' @format A numeric matrix with 10 rows (pathways) and 30 columns (samples).
#' @source Simulated synthetic data for MiCARA testing and demonstration.
"example_pathways"

#' Synthetic Example Sample Metadata
#'
#' A simulated metadata data frame containing sample covariates and disease status.
#'
#' @format A data frame with 30 rows and 6 columns:
#' \describe{
#'   \item{sample_id}{Unique sample identifier}
#'   \item{disease}{Disease group assignment ("Control" or "Case")}
#'   \item{age_category}{Age cohort ("Adult" or "Senior")}
#'   \item{BMI}{Body Mass Index}
#'   \item{gender}{Sex of the participant}
#'   \item{study_name}{Cohort study origin}
#' }
#' @source Simulated synthetic data for MiCARA testing and demonstration.
"example_metadata"
