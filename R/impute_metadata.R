#' Impute Missing Metadata Covariates Using MICE
#'
#' Imputes missing values in host confounders (\code{age}, \code{gender}, \code{BMI})
#' that were retained by \code{\link{validate_inputs}} (i.e., below the missingness
#' cutoff). Experimental design indicators such as \code{disease}, \code{study_name},
#' and sequencing depth (\code{number_reads}) are strictly protected and used only
#' as auxiliary predictors to inform imputation. \code{age_category} is also treated
#' as an auxiliary predictor when complete to guide continuous \code{age} imputation.
#'
#' @details
#' \strong{Scientific Guidance on Downstream Covariate Selection:}
#' Imputed continuous \code{age} is generally recommended over binned \code{age_category}
#' for differential abundance modeling (e.g., ANCOM-BC2) as it avoids discretization
#' power loss. However, \code{age_category} should be used directly if modeling
#' distinct non-linear biological thresholds (e.g., pediatric vs. adult states) or
#' if continuous \code{age} missingness exceeded 50\%.
#'
#' @param micara_obj An object of class \code{micara_input}, as returned
#'   by \code{\link{validate_inputs}}.
#' @param impute_vars Character vector of host confounders to impute.
#'   Defaults to the intersection of \code{c("age", "gender", "BMI")} and
#'   \code{micara_obj$confounders} (i.e., host covariates that passed the
#'   missingness cutoff in \code{validate_inputs()}).
#' @param auxiliary_predictors Character vector of additional metadata
#'   columns to include as MICE predictors (never imputation targets).
#'   Defaults to \code{c("study_name", "age_category", "disease")} when present.
#' @param m Number of multiple imputations (default 5, the MICE convention).
#' @param method Imputation method passed to \code{mice::mice()}. Defaults to
#'   \code{NULL}, allowing MICE to automatically select the optimal algorithm
#'   based on column type (\code{"pmm"} for continuous, \code{"logreg"} for binary,
#'   \code{"polyreg"} for categorical).
#' @param seed Random seed for reproducibility (default 123).
#' @param imputation_index Which of the \code{m} imputed datasets to extract
#'   for the returned \code{metadata} (default 1). The full \code{mids} fit object
#'   is retained in \code{micara_obj$mice_fit} for diagnostic checks or pooled analyses.
#' @param verbose Logical; print progress and summary report to console (default TRUE).
#'
#' @return An object of class \code{c("micara_imputed", "micara_input")}:
#'   the input object with updated \code{metadata}, plus \code{mice_fit} (the full
#'   \code{mids} object, or \code{NULL} if skipped), \code{imputed_vars}, and
#'   \code{imputation_performed} (logical; \code{FALSE} if no missingness was
#'   present or imputation was skipped).
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' sample_ids <- paste0("Sample_", 1:20)
#'
#' taxa_mat <- matrix(
#'     runif(100, min = 0, max = 1000),
#'     nrow = 5, ncol = 20,
#'     dimnames = list(paste0("Taxon_", 1:5), sample_ids)
#' )
#'
#' path_mat <- matrix(
#'     runif(100, min = 0, max = 500),
#'     nrow = 5, ncol = 20,
#'     dimnames = list(paste0("Pathway_", 1:5), sample_ids)
#' )
#'
#' # Create metadata with missing age values (NA)
#' meta_df <- data.frame(
#'     sample = sample_ids,
#'     disease = rep(c("Control", "Disease"), each = 10),
#'     age = c(
#'         25, NA, 30, 45, NA, 50, 33, 29, 41, 38,
#'         52, 60, NA, 48, 55, 62, NA, 39, 44, 51
#'     ),
#'     row.names = sample_ids
#' )
#'
#' # Build micara S3 object
#' validated_obj <- validate_inputs(
#'     taxa = taxa_mat,
#'     pathways = path_mat,
#'     metadata = meta_df,
#'     disease_col = "disease",
#'     disease_min_samples = 5,
#'     verbose = FALSE
#' )
#'
#' # Run imputation on missing age column
#' imputed_obj <- impute_metadata(
#'     micara_obj = validated_obj,
#'     impute_vars = "age",
#'     m = 2,
#'     seed = 123,
#'     verbose = FALSE
#' )
impute_metadata <- function(
  micara_obj,
  impute_vars = NULL,
  auxiliary_predictors = NULL,
  m = 5,
  method = NULL, # (allowing default algorithm selection per data type: pmm (Predictive Mean Matching) for continuous, logreg (Logistic Regression) for binary factors, polyreg (polytomous logistic regression) for unordered factors) or a named vector.
  seed = 123,
  imputation_index = 1,
  verbose = TRUE
) {
    # simple class check
    if (!inherits(micara_obj, "micara_input")) {
        stop(
            "'micara_obj' must be the output of validate_inputs() ",
            "(an object of class 'micara_input')."
        )
    }

    metadata <- micara_obj$metadata

    ## ------------------------------------------------------------------
    ## 1. Determine which variables to impute
    ##
    ## Only true host confounders are ever imputation targets.
    ## age_category is deliberately excluded here - it is expected to be
    ## complete and instead serves as a stable auxiliary predictor for
    ## imputing continuous age (see step 2).
    ## ------------------------------------------------------------------

    host_confounders <- c("age", "gender", "BMI")

    if (is.null(impute_vars)) {
        # Auto-detect intersection of host confounders and retained covariates
        impute_vars <- intersect(host_confounders, micara_obj$confounders)
        # Retain only those that actually contain missing values
        impute_vars <- impute_vars[vapply(impute_vars, function(v) any(is.na(metadata[[v]])), logical(1))]
    } else {
        unavailable <- setdiff(impute_vars, names(metadata))
        if (length(unavailable) > 0) {
            stop("impute_vars not found in metadata: ", paste(unavailable, collapse = ", "))
        }
        if ("age_category" %in% impute_vars) {
            stop(
                "'age_category' should not be used as an imputation target -it is ",
                "expected to be complete and is used as an auxiliary predictor ",
                "instead. Pass it via 'auxiliary_predictors' if you need to include it."
            )
        }
    }
    # Early exit if no variables need imputation
    if (length(impute_vars) == 0) {
        if (verbose) {
            message(
                "[Info] No eligible host confounders with missing values to impute. ",
                "Returning input object unchanged."
            )
        }
        micara_obj$imputation_performed <- FALSE
        micara_obj$mice_fit <- NULL
        micara_obj$imputed_vars <- character(0)
        class(micara_obj) <- c("micara_imputed", class(micara_obj))
        return(micara_obj)
    }

    ## ------------------------------------------------------------------
    ## 2. Determine auxiliary predictors (inform imputation, never imputed)
    ## ------------------------------------------------------------------

    if (is.null(auxiliary_predictors)) {
        auxiliary_predictors <- intersect(
            c("study_name", "age_category", "disease"),
            names(metadata)
        )
    }

    if (any(is.na(metadata[["disease"]]))) {
        stop("'disease' contains missing values; every sample must have a disease label.")
    }

    if ("study_name" %in% auxiliary_predictors && any(is.na(metadata[["study_name"]]))) {
        stop(
            "'study_name' contains missing values. Study/cohort identity should ",
            "never be imputed; either supply it for all samples or remove ",
            "'study_name' from auxiliary_predictors."
        )
    }

    ## age_category is never an imputation target, but it is also not
    ## required to be complete. If it is missing entirely from the
    ## metadata, it is simply absent from auxiliary_predictors already
    ## (see the intersect() above). If it IS present but contains missing
    ## values, we do not block imputation of the actual host confounders
    ## on account of it - we drop it from the predictor set for this run
    ## and proceed using the remaining complete predictors instead.
    if ("age_category" %in% auxiliary_predictors && any(is.na(metadata[["age_category"]]))) {
        if (verbose) {
            message(
                "'age_category' contains missing values and will be excluded as an ",
                "auxiliary predictor for this run (it is never an imputation ",
                "target itself). Imputation of ", paste(impute_vars, collapse = ", "),
                " will proceed using the remaining predictor(s): ",
                paste(setdiff(auxiliary_predictors, "age_category"), collapse = ", "), "."
            )
        }
        auxiliary_predictors <- setdiff(auxiliary_predictors, "age_category")
    }

    ## ------------------------------------------------------------------
    ## 3. Skip if nothing actually needs imputing
    ## ------------------------------------------------------------------

    n_missing <- sum(vapply(metadata[impute_vars], function(x) sum(is.na(x)), integer(1)))

    if (n_missing == 0) {
        if (verbose) {
            message(
                "No missing values detected in ", paste(impute_vars, collapse = ", "),
                "; imputation skipped."
            )
        }
        micara_obj$imputation_performed <- FALSE
        micara_obj$mice_fit <- NULL
        micara_obj$imputed_vars <- impute_vars
        class(micara_obj) <- c("micara_imputed", class(micara_obj))
        return(micara_obj)
    }

    ## ------------------------------------------------------------------
    ## 4. Build the MICE input frame
    ## ------------------------------------------------------------------

    if (!requireNamespace("mice", quietly = TRUE)) {
        stop("Package 'mice' is required for impute_metadata() but is not installed.")
    }

    mice_vars <- unique(c(impute_vars, auxiliary_predictors))
    mice_df <- metadata[, mice_vars, drop = FALSE]

    ## coerce character columns to factors so mice treats them categorically
    mice_df[] <- lapply(mice_df, function(col) {
        if (is.character(col)) factor(col) else col
    })

    ## ------------------------------------------------------------------
    ## 5. Run MICE
    ## ------------------------------------------------------------------

    if (verbose) {
        message("\n[MiCARA] Running MICE imputation for: ", paste(impute_vars, collapse = ", "))
    }

    suppressWarnings({
        mice_fit <- mice::mice(
            mice_df,
            m         = m,
            method    = method,
            seed      = seed,
            printFlag = FALSE
        )
    })

    imputed_df <- mice::complete(mice_fit, imputation_index)

    ## ------------------------------------------------------------------
    ## 6. Merge imputed values back into full metadata
    ## ------------------------------------------------------------------

    for (v in impute_vars) {
        if (is.character(metadata[[v]])) {
            metadata[[v]] <- as.character(imputed_df[[v]])
        } else {
            metadata[[v]] <- imputed_df[[v]]
        }
    }

    metadata[impute_vars] <- imputed_df[impute_vars]

    micara_obj$metadata <- metadata
    micara_obj$mice_fit <- mice_fit
    micara_obj$imputed_vars <- impute_vars
    micara_obj$imputation_performed <- TRUE

    class(micara_obj) <- c("micara_imputed", class(micara_obj))

    ## ------------------------------------------------------------------
    ## 7. Report
    ## ------------------------------------------------------------------

    if (verbose) {
        n_logged <- if (!is.null(mice_fit$loggedEvents)) nrow(mice_fit$loggedEvents) else 0

        imp_summary <- vapply(impute_vars, function(v) {
            n_imp <- mice_fit$nmis[[v]]
            pct_imp <- (n_imp / nrow(metadata)) * 100
            sprintf("  • %-12s : %d / %d samples imputed (%.1f%%)", v, n_imp, nrow(metadata), pct_imp)
        }, FUN.VALUE = character(1))

        logged_text <- if (n_logged > 0) " (see micara_obj$mice_fit$loggedEvents)" else ""

        report <- c(
            "\n========== MiCARA Imputation Report ==========",
            "Imputation Breakdown:",
            imp_summary,
            "------------------------------------------------",
            paste0("Auxiliary predictors:   ", paste(auxiliary_predictors, collapse = ", ")),
            paste0("Imputations (m):        ", m),
            paste0("Dataset used:           #", imputation_index, " of ", m),
            paste0("Logged events:          ", n_logged, logged_text),
            "[MiCARA Tip] Downstream Covariate Selection:",
            "  Imputed continuous 'age' is recommended for linear modeling (ANCOM-BC2)",
            "   rather than 'age_category' to preserve degrees of freedom.",
            "  Use 'age_category' instead if modeling non-linear life-stage thresholds",
            "    (e.g., pediatric vs. adult) or if continuous 'age' missingness was high.",
            "------------------------------------------------"
        )

        message(paste(report, collapse = "\n"))
    }
    return(micara_obj)
}
