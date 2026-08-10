#' Compute confounder-adjusted Residual CLR taxon-pathway correlations for Disease-Associated Features 
#'
#' Applies the second analytical layer of MiCARA after differential
#' abundance analysis. For each disease with at least one robust
#' differentially abundant taxon and pathway, the function:
#'
#' \enumerate{
#'   \item CLR-transforms the complete taxonomic and pathway matrices.
#'   \item selects the disease-specific significant taxa and pathways.
#'   \item restricts the analysis to samples from that disease.
#'   \item residualises both feature sets against estimable covariates.
#'   \item computes all taxon-pathway Spearman correlations on the residuals.
#'   \item applies Benjamini-Hochberg correction within each disease.
#'   \item retains links satisfying both \code{q < q_cutoff} and
#'   \code{abs(rho) > r_cutoff}.
#' }
#'
#' CLR is calculated on the full feature set before disease-specific feature
#' selection. Canonical taxon and pathway identifiers are preserved; labels
#' should be cleaned only in downstream plotting functions.
#'
#' @param micara_obj A \code{micara_input} or \code{micara_imputed} object.
#' @param diffab_obj Optional combined \code{micara_diffab_list} object with
#'   elements \code{$taxa} and \code{$pathways}. A single
#'   \code{micara_diffab} object is also accepted when the complementary
#'   object is supplied through \code{taxa_diffab} or
#'   \code{pathway_diffab}.
#' @param taxa_diffab Optional \code{micara_diffab} object produced by
#'   \code{\link{run_diff_abundance}} with \code{feature_type = "taxa"}.
#' @param pathway_diffab Optional \code{micara_diffab} object produced by
#'   \code{\link{run_diff_abundance}} with
#'   \code{feature_type = "pathways"}.
#' @param diseases Character vector of diseases to analyse. The default
#'   \code{NULL} analyses every disease represented in both differential
#'   abundance result objects.
#' @param covariates Character vector of metadata columns to remove before
#'   correlation. The default \code{NULL} uses retained MiCARA confounders
#'   among \code{age_category}, \code{BMI}, \code{gender}, and
#'   \code{study_name}. This reproduces the CAMDA analysis when all four are
#'   available. Non-estimable covariates are dropped within each disease and
#'   reported.
#' @param q_cutoff Strict BH-adjusted p-value threshold for retaining links.
#'   Default \code{0.05}.
#' @param r_cutoff Strict minimum absolute Spearman correlation for retaining
#'   links. Default \code{0.2}. This is an effect-size threshold for network
#'   inclusion, not a universal biological cutoff.
#' @param pseudocount Zero replacement used before taking logarithms.
#'   \code{NULL} uses half the minimum positive value separately for taxa
#'   and pathways. A single positive number applies to both matrices. A named
#'   list such as \code{list(taxa = 1e-06, pathways = 1e-08)} may also be
#'   supplied.
#' @param min_samples_per_disease Minimum number of disease samples required
#'   before attempting residualisation. Default \code{10}.
#' @param min_residual_df Minimum residual degrees of freedom required after
#'   building the covariate design matrix. Default \code{3}.
#' @param retain_all_pairs Logical. If \code{TRUE}, each disease result also
#'   stores every tested pair and its p/q values in \code{$all_links}.
#'   Default \code{TRUE}.
#' @param verbose Logical; print progress and a disease-level summary.
#'
#' @return An object of class \code{micara_interactions}. It is a named list
#'   with one element per analysed disease. Each disease element contains:
#'   \describe{
#'     \item{status}{\code{"analysed"} or \code{"skipped"}.}
#'     \item{reason}{Reason for skipping, otherwise \code{NA}.}
#'     \item{links}{Associations passing both thresholds.}
#'     \item{all_links}{All tested associations when
#'       \code{retain_all_pairs = TRUE}.}
#'     \item{node_directions}{Differential-abundance directions for tested
#'       taxa and pathways.}
#'     \item{covariates_requested, covariates_used, covariates_dropped}{
#'       Covariate provenance for that disease.}
#'     \item{n_samples, residual_df, n_taxa_tested, n_pathways_tested,
#'       n_pairs_tested}{Analysis scope and model diagnostics.}
#'   }
#'
#'   The object carries \code{"summary"} and \code{"settings"} attributes.
#'   Access them with \code{attr(x, "summary")} and
#'   \code{attr(x, "settings")}.
#'
#' @export
compute_residualCLR_correlations <- function(
    micara_obj,
    diffab_obj              = NULL,
    taxa_diffab             = NULL,
    pathway_diffab          = NULL,
    diseases                = NULL,
    covariates              = NULL,
    q_cutoff                = 0.05,
    r_cutoff                = 0.2,
    pseudocount             = NULL,
    min_samples_per_disease = 10L,
    min_residual_df         = 3L,
    retain_all_pairs        = TRUE,
    verbose                 = TRUE
) {
  

  ## ------------------------------------------------------------------
  ## 1. Validate inputs and resolve differential-abundance objects
  ## ------------------------------------------------------------------
  
  #check class types
  if (!inherits(micara_obj, "micara_input")) {
    stop(
      "'micara_obj' must inherit from 'micara_input'. Run ",
      "validate_inputs() first and impute_metadata() when required.",
      call. = FALSE
    )
  }
  # Hmisc availability
  if (!requireNamespace("Hmisc", quietly = TRUE)) {
    stop(
      "Package 'Hmisc' is required. Install it before running ",
      "compute_residual_correlations().",
      call. = FALSE
    )
  }
  #Validates that thresholds fall within acceptable mathematical ranges
  .check_probability(q_cutoff, "q_cutoff", allow_zero = FALSE)
  .check_probability(r_cutoff, "r_cutoff", allow_zero = TRUE)

  if (!is.numeric(min_samples_per_disease) ||
      length(min_samples_per_disease) != 1L ||
      is.na(min_samples_per_disease) ||
      min_samples_per_disease < 3) {
    stop("'min_samples_per_disease' must be one number >= 3.", call. = FALSE)
  }

  if (!is.numeric(min_residual_df) ||
      length(min_residual_df) != 1L ||
      is.na(min_residual_df) ||
      min_residual_df < 1) {
    stop("'min_residual_df' must be one positive number.", call. = FALSE)
  }

  if (!is.logical(retain_all_pairs) || length(retain_all_pairs) != 1L ||
      is.na(retain_all_pairs)) {
    stop("'retain_all_pairs' must be TRUE or FALSE.", call. = FALSE)
  }
  ## ----------------------------------------------------------------
  ## 2. Metadata & differential abundance resolution
  ## ----------------------------------------------------------------
  
  resolved <- .resolve_diffab_objects(
    diffab_obj     = diffab_obj,
    taxa_diffab    = taxa_diffab,
    pathway_diffab = pathway_diffab
  )
  taxa_diffab    <- resolved$taxa
  pathway_diffab <- resolved$pathways

  .validate_diffab_object(taxa_diffab, "taxa")
  .validate_diffab_object(pathway_diffab, "pathways")

  metadata <- as.data.frame(micara_obj$metadata, stringsAsFactors = FALSE)
  disease_col   <- attr(micara_obj, "disease_col") %||% "disease"
  sample_id_col <- attr(micara_obj, "sample_id_col") %||% "sample_id"
  
  #Checks Metadata
  required_metadata <- c("sample_id", "disease")
  missing_metadata <- setdiff(required_metadata, names(metadata))
  if (length(missing_metadata) > 0L) {
    stop(
      "Missing required metadata column(s): ",
      paste(missing_metadata, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyNA(metadata$sample_id) || any(trimws(as.character(metadata$sample_id)) == "")) {
    stop("'metadata$sample_id' contains missing or empty values.", call. = FALSE)
  }
  if (anyDuplicated(metadata$sample_id)) {
    stop("'metadata$sample_id' must contain unique sample IDs.", call. = FALSE)
  }
  if (anyNA(metadata$disease) || any(trimws(as.character(metadata$disease)) == "")) {
    stop("'metadata$disease' contains missing or empty values.", call. = FALSE)
  }

  metadata$sample_id <- as.character(metadata$sample_id)
  metadata$disease   <- make.names(as.character(metadata$disease))

  ## ------------------------------------------------------------------
  ## 3. Resolve & Validate Covariates
  ## ------------------------------------------------------------------
  if (is.null(covariates)) {
    preferred <- c("age_category", "BMI", "gender", "study_name")
    retained  <- micara_obj$confounders %||% character(0)
    covariates <- intersect(preferred, retained)
    
    # Dynamic fallback: grab non-ID/non-disease columns if retain list is empty
    if (length(covariates) == 0L) {
      exclude_cols <- c(sample_id_col, disease_col, "study_name")
      covariates   <- setdiff(colnames(metadata), exclude_cols)
    }
    
    if (verbose && length(covariates) > 0L) {
      message("[Info] Auto-selected covariates for residualization: ", paste(covariates, collapse = ", "))
    }
  }
  
  covariates <- unique(as.character(covariates))
  covariates <- setdiff(covariates, c(sample_id_col, disease_col))
  
  missing_covariates <- setdiff(covariates, names(metadata))
  if (length(missing_covariates) > 0L) {
    stop("Covariate(s) not found in metadata: ", paste(missing_covariates, collapse = ", "), call. = FALSE)
  }
  
  # Validate missingness on covariates
  if (length(covariates) > 0L) {
    cov_df <- metadata[, covariates, drop = FALSE]
    if (anyNA(cov_df)) {
      na_cols <- covariates[colSums(is.na(cov_df)) > 0]
      stop(sprintf(
        "Missing values (NA) in covariates: %s. Run 'impute_metadata()' before running this step.",
        paste(na_cols, collapse = ", ")
      ), call. = FALSE)
    }
  }

  ## ------------------------------------------------------------------
  ## 4. CLR-transform the complete feature matrices once
  ## ------------------------------------------------------------------

  pseudocounts <- .resolve_pseudocounts(pseudocount)

  if (verbose) {
    message("[MiCARA] CLR-transforming complete taxa and pathway matrices.")
  }

  taxa_clr_result <- .clr_transform(
    micara_obj$taxa,
    pseudocount = pseudocounts$taxa,
    label = "taxa"
  )
  pathway_clr_result <- .clr_transform(
    micara_obj$pathways,
    pseudocount = pseudocounts$pathways,
    label = "pathways"
  )

  clr_taxa     <- taxa_clr_result$matrix
  clr_pathways <- pathway_clr_result$matrix

  ## ------------------------------------------------------------------
  ## 3. Determine disease scope
  ## ------------------------------------------------------------------

  taxa_diseases <- unique(as.character(taxa_diffab$significant$disease))
  path_diseases <- unique(as.character(pathway_diffab$significant$disease))
  taxa_diseases <- taxa_diseases[!is.na(taxa_diseases)]
  path_diseases <- path_diseases[!is.na(path_diseases)]

  diseases_with_both <- intersect(taxa_diseases, path_diseases)
  skipped_no_pathways <- setdiff(taxa_diseases, path_diseases)
  skipped_no_taxa     <- setdiff(path_diseases, taxa_diseases)

  if (is.null(diseases)) {
    diseases <- diseases_with_both
  } else {
    diseases <- unique(as.character(diseases))

    unavailable <- setdiff(diseases, diseases_with_both)
    if (length(unavailable) > 0L) {
      warning(
        "The following requested disease(s) do not have both significant ",
        "taxa and significant pathways and will not be analysed: ",
        paste(unavailable, collapse = ", "),
        call. = FALSE
      )
    }
    diseases <- intersect(diseases, diseases_with_both)
  }

  if (verbose && length(skipped_no_pathways) > 0L) {
    message(
      "[MiCARA] Significant taxa but no significant pathways: ",
      paste(skipped_no_pathways, collapse = ", "),
      ". These diseases are not eligible for interaction inference."
    )
  }
  if (verbose && length(skipped_no_taxa) > 0L) {
    message(
      "[MiCARA] Significant pathways but no significant taxa: ",
      paste(skipped_no_taxa, collapse = ", "),
      ". These diseases are not eligible for interaction inference."
    )
  }

  if (length(diseases) == 0L) {
    stop(
      "No diseases have both significant taxa and significant pathways.",
      call. = FALSE
    )
  }

  ## ------------------------------------------------------------------
  ## 4. Run each disease independently
  ## ------------------------------------------------------------------

  results <- lapply(diseases, function(disease_name) {
    .process_one_disease(
      disease                  = disease_name,
      metadata                 = metadata,
      clr_taxa                 = clr_taxa,
      clr_pathways             = clr_pathways,
      taxa_diffab              = taxa_diffab,
      pathway_diffab           = pathway_diffab,
      covariates               = covariates,
      q_cutoff                 = q_cutoff,
      r_cutoff                 = r_cutoff,
      min_samples_per_disease  = as.integer(min_samples_per_disease),
      min_residual_df          = as.integer(min_residual_df),
      retain_all_pairs         = retain_all_pairs,
      verbose                  = verbose
    )
  })
  names(results) <- diseases

  summary_table <- do.call(
    rbind,
    lapply(names(results), function(disease_name) {
      x <- results[[disease_name]]
      data.frame(
        disease              = disease_name,
        status               = x$status,
        reason               = if (is.null(x$reason)) NA_character_ else x$reason,
        n_samples            = x$n_samples,
        residual_df          = x$residual_df,
        n_taxa_tested        = x$n_taxa_tested,
        n_pathways_tested    = x$n_pathways_tested,
        n_pairs_tested       = x$n_pairs_tested,
        n_significant_links  = nrow(x$links),
        n_positive_links     = sum(x$links$direction == "positive"),
        n_negative_links     = sum(x$links$direction == "negative"),
        covariates_used      = paste(x$covariates_used, collapse = ", "),
        covariates_dropped   = paste(x$covariates_dropped, collapse = ", "),
        stringsAsFactors     = FALSE
      )
    })
  )
  rownames(summary_table) <- NULL

  settings <- list(
    covariates_requested    = covariates,
    cor_method              = "spearman",
    p_adjust_method         = "BH_within_disease",
    q_cutoff                = q_cutoff,
    r_cutoff                = r_cutoff,
    threshold_rule          = "q < q_cutoff and abs(rho) > r_cutoff",
    taxa_pseudocount        = taxa_clr_result$pseudocount,
    pathway_pseudocount     = pathway_clr_result$pseudocount,
    taxa_all_zero_removed   = taxa_clr_result$removed_all_zero_features,
    pathway_all_zero_removed = pathway_clr_result$removed_all_zero_features,
    min_samples_per_disease = as.integer(min_samples_per_disease),
    min_residual_df         = as.integer(min_residual_df),
    retain_all_pairs        = retain_all_pairs,
    skipped_no_pathways     = skipped_no_pathways,
    skipped_no_taxa         = skipped_no_taxa
  )

  attr(results, "summary")  <- summary_table
  attr(results, "settings") <- settings
  class(results) <- c("micara_interactions", "list")

  if (verbose) {
    cat("\n========== MiCARA Residual Correlation Summary ==========\n")
    print(summary_table, row.names = FALSE)
    cat("----------------------------------------------------------\n")
  }

  results
}


## -------------------------------------------------------------------------
## Internal helpers
## -------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.check_probability <- function(x, name, allow_zero = TRUE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    lower_text <- if (allow_zero) "[0, 1]" else "(0, 1]"
    stop("'", name, "' must be a single number in ", lower_text, ".", call. = FALSE)
  }

  lower_ok <- if (allow_zero) x >= 0 else x > 0

  if (!lower_ok || x > 1) {
    lower_text <- if (allow_zero) "[0, 1]" else "(0, 1]"
    stop("'", name, "' must be a single number in ", lower_text, ".", call. = FALSE)
  }

  invisible(TRUE)
}


#' @keywords internal
#' @noRd
.resolve_diffab_objects <- function(diffab_obj, taxa_diffab, pathway_diffab) {

  if (!is.null(diffab_obj) && inherits(diffab_obj, "micara_diffab_list")) {
    if (!is.null(taxa_diffab) || !is.null(pathway_diffab)) {
      stop(
        "Supply either a combined 'diffab_obj' or separate taxa/pathway ",
        "objects, not both.",
        call. = FALSE
      )
    }

    taxa_diffab    <- diffab_obj$taxa
    pathway_diffab <- diffab_obj$pathways
  } else if (!is.null(diffab_obj) && inherits(diffab_obj, "micara_diffab")) {
    if (identical(diffab_obj$feature_type, "taxa")) {
      if (!is.null(taxa_diffab)) {
        stop("Taxa differential-abundance results were supplied twice.", call. = FALSE)
      }
      taxa_diffab <- diffab_obj
    } else if (identical(diffab_obj$feature_type, "pathways")) {
      if (!is.null(pathway_diffab)) {
        stop("Pathway differential-abundance results were supplied twice.", call. = FALSE)
      }
      pathway_diffab <- diffab_obj
    } else {
      stop("'diffab_obj$feature_type' must be 'taxa' or 'pathways'.", call. = FALSE)
    }
  } else if (!is.null(diffab_obj)) {
    stop(
      "'diffab_obj' must be a 'micara_diffab' or 'micara_diffab_list' object.",
      call. = FALSE
    )
  }

  list(taxa = taxa_diffab, pathways = pathway_diffab)
}


#' @keywords internal
#' @noRd
.validate_diffab_object <- function(x, expected_type) {

  if (is.null(x) || !inherits(x, "micara_diffab")) {
    stop(
      "A valid 'micara_diffab' object is required for ", expected_type, ".",
      call. = FALSE
    )
  }

  if (!identical(x$feature_type, expected_type)) {
    stop(
      "Expected feature_type = '", expected_type, "' but received '",
      x$feature_type, "'.",
      call. = FALSE
    )
  }

  required_columns <- c("feature", "disease", "lfc")
  missing_columns <- setdiff(required_columns, names(x$significant))
  if (length(missing_columns) > 0L) {
    stop(
      "The ", expected_type, " differential-abundance object's ",
      "'significant' table is missing: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' @keywords internal
#' @noRd
.resolve_pseudocounts <- function(pseudocount) {

  if (is.null(pseudocount)) {
    return(list(taxa = NULL, pathways = NULL))
  }

  if (is.numeric(pseudocount) && length(pseudocount) == 1L &&
      is.finite(pseudocount) && pseudocount > 0) {
    return(list(taxa = pseudocount, pathways = pseudocount))
  }

  if (is.list(pseudocount) &&
      all(c("taxa", "pathways") %in% names(pseudocount))) {

    valid_value <- function(x) {
      is.null(x) ||
        (is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0)
    }

    if (!valid_value(pseudocount$taxa) ||
        !valid_value(pseudocount$pathways)) {
      stop(
        "Named pseudocount values must be NULL or single positive numbers.",
        call. = FALSE
      )
    }

    return(list(
      taxa     = pseudocount$taxa,
      pathways = pseudocount$pathways
    ))
  }

  stop(
    "'pseudocount' must be NULL, one positive number, or a named list ",
    "with elements 'taxa' and 'pathways'.",
    call. = FALSE
  )
}


#' @keywords internal
#' @noRd
.clr_transform <- function(mat, pseudocount = NULL, label = "matrix") {

  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"

  if (is.null(rownames(mat)) || any(rownames(mat) == "")) {stop("The ", label, " matrix must have non-empty feature names.", call. = FALSE)}
  if (is.null(colnames(mat)) || any(colnames(mat) == "")) {stop("The ", label, " matrix must have non-empty sample IDs.", call. = FALSE)}
  if (anyDuplicated(rownames(mat))) {stop("The ", label, " matrix contains duplicated feature names.", call. = FALSE)}
  if (anyDuplicated(colnames(mat))) {stop("The ", label, " matrix contains duplicated sample IDs.", call. = FALSE)}
  if (anyNA(mat) || any(!is.finite(mat))) {stop("The ", label, " matrix contains missing or non-finite values.", call. = FALSE)}
  if (any(mat < 0)) {stop("The ", label, " matrix contains negative values and cannot be CLR-transformed.", call. = FALSE)
  }

  all_zero_features <- rownames(mat)[rowSums(mat > 0) == 0L]
  if (length(all_zero_features) > 0L) {
    warning(
      "Removing ", length(all_zero_features), " all-zero ", label,
      " feature(s) before CLR transformation.",
      call. = FALSE
    )
    mat <- mat[rowSums(mat > 0) > 0L, , drop = FALSE]
  }

  if (nrow(mat) == 0L) {
    stop(
      "The ", label, " matrix contains no informative features after ",
      "removing all-zero rows.",
      call. = FALSE
    )
  }

  all_zero_samples <- colnames(mat)[colSums(mat > 0) == 0L]
  if (length(all_zero_samples) > 0L) {
    warning(
      "Removing ", length(all_zero_samples), " all-zero sample(s) from ", label,
      " matrix: ", paste(utils::head(all_zero_samples, 5L), collapse = ", "),
      if (length(all_zero_samples) > 5L) " ..." else "",
      call. = FALSE
    )
    mat <- mat[, colSums(mat > 0) > 0L, drop = FALSE]
  }
  
  if (ncol(mat) == 0L) {
    stop(
      "The ", label, " matrix contains no samples after removing all-zero columns.",
      call. = FALSE
    )
  }
  
  # 3. Calculate pseudocount & CLR
  if (is.null(pseudocount)) {
    positive_values <- mat[mat > 0]
    if (length(positive_values) == 0L) {
      stop("The ", label, " matrix has no positive values.", call. = FALSE)
    }
    pseudocount <- min(positive_values) / 2
  }

  if (!is.numeric(pseudocount) || length(pseudocount) != 1L ||
      is.na(pseudocount) || !is.finite(pseudocount) || pseudocount <= 0) {
    stop("The ", label, " pseudocount must be one positive number.", call. = FALSE)
  }

  mat[mat == 0] <- pseudocount
  log_mat <- log(mat)
  clr_mat <- sweep(log_mat, 2L, colMeans(log_mat), FUN = "-")

  list(
    matrix = clr_mat,
    pseudocount = pseudocount,
    removed_all_zero_features = all_zero_features
  )
}


#' @keywords internal
#' @noRd
.prepare_residual_design <- function(
    metadata,
    covariates,
    min_residual_df
) {

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  covariates_requested <- covariates
  covariates_used <- character(0)
  covariates_dropped <- character(0)

  for (covariate in covariates_requested) {
    values <- metadata[[covariate]]

    if (anyNA(values)) {
      return(list(
        ok = FALSE,
        reason = paste0("missing values in covariate '", covariate, "'"),
        covariates_used = covariates_used,
        covariates_dropped = covariates_dropped,
        residual_df = NA_integer_
      ))
    }

    if (is.character(values) || is.logical(values)) {
      values <- factor(values)
    }
    if (is.factor(values)) {
      values <- droplevels(values)
    }

    metadata[[covariate]] <- values

    if (length(unique(values)) <= 1L) {
      covariates_dropped <- c(
        covariates_dropped,
        paste0(covariate, " [zero variance]")
      )
    } else {
      covariates_used <- c(covariates_used, covariate)
    }
  }

  formula_object <- if (length(covariates_used) == 0L) {
    stats::as.formula("~ 1")
  } else {
    stats::reformulate(covariates_used)
  }

  design <- stats::model.matrix(formula_object, data = metadata)
  design_qr <- qr(design)
  design_rank <- design_qr$rank
  residual_df <- nrow(design) - design_rank

  aliased_columns <- character(0)
  if (design_rank < ncol(design)) {
    independent <- design_qr$pivot[seq_len(design_rank)]
    aliased <- setdiff(seq_len(ncol(design)), independent)
    aliased_columns <- colnames(design)[aliased]
    covariates_dropped <- c(
      covariates_dropped,
      paste0(aliased_columns, " [aliased design column]")
    )
  }

  if (residual_df < min_residual_df) {
    return(list(
      ok = FALSE,
      reason = paste0(
        "only ", residual_df, " residual degree(s) of freedom; minimum is ",
        min_residual_df
      ),
      covariates_used = covariates_used,
      covariates_dropped = covariates_dropped,
      residual_df = residual_df,
      design_rank = design_rank,
      aliased_columns = aliased_columns
    ))
  }

  list(
    ok = TRUE,
    reason = NA_character_,
    metadata = metadata,
    design = design,
    design_qr = design_qr,
    design_rank = design_rank,
    residual_df = residual_df,
    covariates_used = covariates_used,
    covariates_dropped = unique(covariates_dropped),
    aliased_columns = aliased_columns
  )
}


#' @keywords internal
#' @noRd
.residualize_with_design <- function(mat, design_qr) {

  mat <- as.matrix(mat)
  residuals <- qr.resid(design_qr, t(mat))
  residuals <- t(residuals)

  rownames(residuals) <- rownames(mat)
  colnames(residuals) <- colnames(mat)

  residuals
}


#' @keywords internal
#' @noRd
.remove_constant_rows <- function(mat, tolerance = sqrt(.Machine$double.eps)) {

  row_variances <- apply(mat, 1L, stats::var)
  keep <- is.finite(row_variances) & row_variances > tolerance

  list(
    matrix = mat[keep, , drop = FALSE],
    removed = rownames(mat)[!keep]
  )
}


#' @keywords internal
#' @noRd
.spearman_pairs <- function(taxa_residuals, pathway_residuals) {

  combined <- rbind(taxa_residuals, pathway_residuals)
  correlation_result <- Hmisc::rcorr(t(combined), type = "spearman")

  n_taxa <- nrow(taxa_residuals)
  n_pathways <- nrow(pathway_residuals)

  pathway_indices <- n_taxa + seq_len(n_pathways)

  list(
    rho = correlation_result$r[
      seq_len(n_taxa),
      pathway_indices,
      drop = FALSE
    ],
    p = correlation_result$P[
      seq_len(n_taxa),
      pathway_indices,
      drop = FALSE
    ]
  )
}


#' @keywords internal
#' @noRd
.empty_links_table <- function() {
  data.frame(
    taxon = character(0),
    pathway = character(0),
    rho = numeric(0),
    p = numeric(0),
    q = numeric(0),
    direction = character(0),
    stringsAsFactors = FALSE
  )
}


#' @keywords internal
#' @noRd
.skipped_disease_result <- function(
    reason,
    covariates_requested,
    covariates_used = character(0),
    covariates_dropped = character(0),
    n_samples = 0L,
    residual_df = NA_integer_,
    n_taxa_tested = 0L,
    n_pathways_tested = 0L
) {
  list(
    status = "skipped",
    reason = reason,
    links = .empty_links_table(),
    all_links = NULL,
    node_directions = data.frame(
      name = character(0),
      type = character(0),
      group = character(0),
      stringsAsFactors = FALSE
    ),
    covariates_requested = covariates_requested,
    covariates_used = covariates_used,
    covariates_dropped = covariates_dropped,
    missing_taxa = character(0),
    missing_pathways = character(0),
    constant_taxa = character(0),
    constant_pathways = character(0),
    n_samples = as.integer(n_samples),
    residual_df = as.integer(residual_df),
    n_taxa_tested = as.integer(n_taxa_tested),
    n_pathways_tested = as.integer(n_pathways_tested),
    n_pairs_tested = 0L
  )
}


#' @keywords internal
#' @noRd
.process_one_disease <- function(
    disease,
    metadata,
    clr_taxa,
    clr_pathways,
    taxa_diffab,
    pathway_diffab,
    covariates,
    q_cutoff,
    r_cutoff,
    min_samples_per_disease,
    min_residual_df,
    retain_all_pairs,
    verbose
) {

  taxa_significant <- taxa_diffab$significant[
    taxa_diffab$significant$disease == disease,
    ,
    drop = FALSE
  ]
  pathway_significant <- pathway_diffab$significant[
    pathway_diffab$significant$disease == disease,
    ,
    drop = FALSE
  ]

  significant_taxa <- unique(as.character(taxa_significant$feature))
  significant_pathways <- unique(as.character(pathway_significant$feature))

  missing_taxa <- setdiff(significant_taxa, rownames(clr_taxa))
  missing_pathways <- setdiff(significant_pathways, rownames(clr_pathways))

  if (length(missing_taxa) > 0L) {
    warning(
      disease, ": ", length(missing_taxa),
      " significant taxon/taxa were absent from the taxa matrix and were dropped.",
      call. = FALSE
    )
  }
  if (length(missing_pathways) > 0L) {
    warning(
      disease, ": ", length(missing_pathways),
      " significant pathway(s) were absent from the pathway matrix and were dropped.",
      call. = FALSE
    )
  }

  significant_taxa <- intersect(significant_taxa, rownames(clr_taxa))
  significant_pathways <- intersect(significant_pathways, rownames(clr_pathways))

  if (length(significant_taxa) == 0L || length(significant_pathways) == 0L) {
    return(.skipped_disease_result(
      reason = "no usable significant taxa or pathways after feature alignment",
      covariates_requested = covariates,
      n_taxa_tested = length(significant_taxa),
      n_pathways_tested = length(significant_pathways)
    ))
  }

  disease_samples <- metadata$sample_id[metadata$disease == disease]
  common_samples <- Reduce(
    intersect,
    list(disease_samples, colnames(clr_taxa), colnames(clr_pathways))
  )

  if (length(common_samples) < min_samples_per_disease) {
    return(.skipped_disease_result(
      reason = paste0(
        "only ", length(common_samples), " aligned sample(s); minimum is ",
        min_samples_per_disease
      ),
      covariates_requested = covariates,
      n_samples = length(common_samples),
      n_taxa_tested = length(significant_taxa),
      n_pathways_tested = length(significant_pathways)
    ))
  }

  metadata_index <- match(common_samples, metadata$sample_id)
  if (anyNA(metadata_index)) {
    return(.skipped_disease_result(
      reason = "sample alignment failed while matching metadata",
      covariates_requested = covariates,
      n_samples = length(common_samples),
      n_taxa_tested = length(significant_taxa),
      n_pathways_tested = length(significant_pathways)
    ))
  }

  metadata_subset <- metadata[metadata_index, , drop = FALSE]
  rownames(metadata_subset) <- common_samples

  if (!identical(metadata_subset$sample_id, common_samples)) {
    return(.skipped_disease_result(
      reason = "metadata sample order does not match abundance matrices",
      covariates_requested = covariates,
      n_samples = length(common_samples),
      n_taxa_tested = length(significant_taxa),
      n_pathways_tested = length(significant_pathways)
    ))
  }

  design_result <- .prepare_residual_design(
    metadata = metadata_subset,
    covariates = covariates,
    min_residual_df = min_residual_df
  )

  if (!design_result$ok) {
    if (verbose) {
      message("[MiCARA] ", disease, " skipped: ", design_result$reason, ".")
    }

    return(.skipped_disease_result(
      reason = design_result$reason,
      covariates_requested = covariates,
      covariates_used = design_result$covariates_used,
      covariates_dropped = design_result$covariates_dropped,
      n_samples = length(common_samples),
      residual_df = design_result$residual_df,
      n_taxa_tested = length(significant_taxa),
      n_pathways_tested = length(significant_pathways)
    ))
  }

  taxa_matrix <- clr_taxa[
    significant_taxa,
    common_samples,
    drop = FALSE
  ]
  pathway_matrix <- clr_pathways[
    significant_pathways,
    common_samples,
    drop = FALSE
  ]

  taxa_residuals <- .residualize_with_design(
    taxa_matrix,
    design_result$design_qr
  )
  pathway_residuals <- .residualize_with_design(
    pathway_matrix,
    design_result$design_qr
  )

  taxa_constant_check <- .remove_constant_rows(taxa_residuals)
  pathway_constant_check <- .remove_constant_rows(pathway_residuals)

  taxa_residuals <- taxa_constant_check$matrix
  pathway_residuals <- pathway_constant_check$matrix

  if (nrow(taxa_residuals) == 0L || nrow(pathway_residuals) == 0L) {
    return(.skipped_disease_result(
      reason = "all residualised taxa or pathways had zero variance",
      covariates_requested = covariates,
      covariates_used = design_result$covariates_used,
      covariates_dropped = design_result$covariates_dropped,
      n_samples = length(common_samples),
      residual_df = design_result$residual_df,
      n_taxa_tested = nrow(taxa_residuals),
      n_pathways_tested = nrow(pathway_residuals)
    ))
  }

  pair_statistics <- .spearman_pairs(
    taxa_residuals,
    pathway_residuals
  )

  all_links <- data.frame(
    taxon = rep(
      rownames(pair_statistics$rho),
      times = ncol(pair_statistics$rho)
    ),
    pathway = rep(
      colnames(pair_statistics$rho),
      each = nrow(pair_statistics$rho)
    ),
    rho = as.vector(pair_statistics$rho),
    p = as.vector(pair_statistics$p),
    stringsAsFactors = FALSE
  )

  all_links$q <- stats::p.adjust(all_links$p, method = "BH")
  all_links$direction <- ifelse(
    is.na(all_links$rho),
    NA_character_,
    ifelse(all_links$rho > 0, "positive",
           ifelse(all_links$rho < 0, "negative", "zero"))
  )

  links <- all_links[
    !is.na(all_links$q) &
      all_links$q < q_cutoff &
      abs(all_links$rho) > r_cutoff,
    ,
    drop = FALSE
  ]
  rownames(links) <- NULL

  taxa_lfc_lookup <- stats::setNames(
    taxa_significant$lfc,
    taxa_significant$feature
  )
  pathway_lfc_lookup <- stats::setNames(
    pathway_significant$lfc,
    pathway_significant$feature
  )

  direction_group <- function(values) {
    ifelse(
      is.na(values),
      "unknown",
      ifelse(values > 0, "up", ifelse(values < 0, "down", "unchanged"))
    )
  }

  node_directions <- rbind(
    data.frame(
      name = rownames(taxa_residuals),
      type = "Taxon",
      group = direction_group(
        taxa_lfc_lookup[rownames(taxa_residuals)]
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      name = rownames(pathway_residuals),
      type = "Pathway",
      group = direction_group(
        pathway_lfc_lookup[rownames(pathway_residuals)]
      ),
      stringsAsFactors = FALSE
    )
  )
  rownames(node_directions) <- NULL

  if (verbose) {
    message(
      "[MiCARA] ", disease, ": ",
      nrow(taxa_residuals), " taxa x ",
      nrow(pathway_residuals), " pathways = ",
      nrow(all_links), " tested pair(s); ",
      nrow(links), " retained link(s)."
    )
  }

  list(
    status = "analysed",
    reason = NA_character_,
    links = links,
    all_links = if (retain_all_pairs) all_links else NULL,
    node_directions = node_directions,
    covariates_requested = covariates,
    covariates_used = design_result$covariates_used,
    covariates_dropped = design_result$covariates_dropped,
    missing_taxa = missing_taxa,
    missing_pathways = missing_pathways,
    constant_taxa = taxa_constant_check$removed,
    constant_pathways = pathway_constant_check$removed,
    n_samples = length(common_samples),
    residual_df = design_result$residual_df,
    n_taxa_tested = nrow(taxa_residuals),
    n_pathways_tested = nrow(pathway_residuals),
    n_pairs_tested = nrow(all_links)
  )
}
