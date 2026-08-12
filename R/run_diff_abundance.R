#' Run Differential Abundance Analysis Using ANCOM-BC2
#'
#' Performs bias-corrected differential abundance testing for microbial taxa,
#' functional pathways, or both using ANCOM-BC2 (\code{ANCOMBC::ancombc2}).
#' Automatically integrates retained host confounders from \code{\link{validate_inputs}}
#' as fixed effects and supports multi-cohort random effects.
#'
#' @details
#' \strong{Multi-Feature Execution:}
#' When \code{feature_type = c("taxa", "pathways")}, the function iterates over
#' both datasets if present in \code{micara_obj} and returns a named list of
#' analysis objects.
#'
#' \strong{Model Specification & Safeguards:}
#' Fixed effects automatically incorporate retained host confounders (e.g., age, gender,
#' BMI) alongside \code{disease}. Variables with zero variance or missing values are
#' automatically flagged. The random effect (\code{rand_effect = "study_name"})
#' accounts for cohort-level batch variation and is automatically disabled if only a
#' single cohort/study is present in the dataset.
#'
#' \strong{Sensitivity Flags:}
#' By default (\code{use_robust = TRUE}), significance is defined using robust sensitivity
#' flags (\code{diff_robust_disease...}) calculated across pseudo-count additions,
#' providing high confidence against false positives in sparse microbiome data.
#'
#' @param micara_obj An object of class \code{micara_input} (or \code{micara_imputed}),
#'   as returned by \code{\link{validate_inputs}} or \code{\link{impute_metadata}}.
#' @param feature_type Character vector specifying feature matrices to analyze:
#'   \code{"taxa"}, \code{"pathways"}, or \code{c("taxa", "pathways")} (default).
#'   Keywords \code{"all"} and \code{"both"} are also supported.
#' @param fix_effects Character vector of metadata columns to include as fixed-effect
#'   covariates. If \code{NULL} (default), auto-detects retained confounders from
#'   \code{micara_obj$confounders} (excluding \code{"study_name"} and \code{"disease"}).
#' @param exclude_covariates Character vector of covariate names to explicitly exclude
#'   from fixed effects.
#' @param rand_effect Character string naming the metadata column for random intercepts
#'   (default \code{"study_name"}). Automatically set to \code{NULL} if the column is
#'   missing or contains only one unique level.
#' @param reference_level Character string defining the control or reference category
#'   within \code{metadata$disease} (default \code{"healthy"}). All differential
#'   log-fold changes will be calculated relative to this baseline.
#' @param p_adj_method Method for p-value adjustment passed to \code{ancombc2}
#'   (default \code{"BH"}).
#' @param prv_cut Minimum prevalence threshold (0 to 1) for features to be retained
#'   prior to analysis (default 0).
#' @param alpha Significance level threshold for q-values (default 0.05).
#' @param pairwise Logical; whether to perform pairwise comparisons between all disease
#'   levels (default FALSE).
#' @param struc_zero Logical; whether to detect structural zeros (default TRUE).
#' @param neg_lb Logical; whether to apply a negative lower bound for asymptotic variance
#'   estimation (default TRUE).
#' @param pseudo_sens Logical; whether to perform pseudo-count sensitivity analysis
#'   (default TRUE). Automatically set to \code{TRUE} if \code{use_robust = TRUE}.
#' @param use_robust Logical; whether to use robust sensitivity flags
#'   (\code{diff_robust_disease}) rather than standard differential abundance flags
#'   (\code{diff_disease}) to define significance (default TRUE).
#' @param n_cores Integer specifying the number of CPU cores for parallel processing
#'   (default 1).
#' @param verbose Logical; print progress and summary report to console (default TRUE).
#'
#' @return If a single \code{feature_type} is analyzed, an object of class
#'   \code{"micara_diffab"}. If multiple feature types are analyzed, a named list
#'   of \code{"micara_diffab"} objects with class \code{"micara_diffab_list"}.
#'
#' @importFrom phyloseq phyloseq otu_table sample_data
#' @importFrom ANCOMBC ancombc2
#' @importFrom stats relevel na.omit
#' @export
#'
#' @examples
#' set.seed(123)
#' sample_ids <- paste0("Sample_", 1:20)
#'
#' # Create count matrices required for differential abundance modeling
#' taxa_mat <- matrix(
#'     rpois(100, lambda = 50),
#'     nrow = 5, ncol = 20,
#'     dimnames = list(paste0("Taxon_", 1:5), sample_ids)
#' )
#'
#' path_mat <- matrix(
#'     rpois(100, lambda = 100),
#'     nrow = 5, ncol = 20,
#'     dimnames = list(paste0("Pathway_", 1:5), sample_ids)
#' )
#'
#' meta_df <- data.frame(
#'     sample = sample_ids,
#'     disease = rep(c("Control", "Disease"), each = 10),
#'     study_name = rep(c("Study1", "Study2"), each = 10),
#'     row.names = sample_ids
#' )
#'
#' validated_obj <- validate_inputs(
#'     taxa = taxa_mat,
#'     pathways = path_mat,
#'     metadata = meta_df,
#'     disease_col = "disease",
#'     disease_min_samples = 5,
#'     verbose = FALSE
#' )
#'
#' # Run differential abundance analysis (set n_cores = 1 for CRAN/Bioc compliance)
#' diff_res <- run_diff_abundance(
#'     micara_obj = validated_obj,
#'     feature_type = "taxa",
#'     rand_effect = NULL,
#'     n_cores = 1,
#'     verbose = FALSE
#' )
run_diff_abundance <- function(
  micara_obj,
  feature_type = c("taxa", "pathways"),
  fix_effects = NULL,
  exclude_covariates = NULL,
  rand_effect = "study_name",
  reference_level = NULL, # Set default to NULL for dynamic matching
  p_adj_method = "BH",
  prv_cut = 0.10, # Set to 10% prevalence to prevent infinite fitting loops
  alpha = 0.05,
  pairwise = FALSE,
  struc_zero = TRUE,
  neg_lb = TRUE,
  pseudo_sens = TRUE,
  use_robust = TRUE,
  n_cores = 1, # Primary argument for cores
  verbose = TRUE
) {
    # Auto-detection of reference level terms without throwing an error
    # Detect disease column and unique levels
    disease_col <- if (!is.null(attr(micara_obj, "disease_col"))) attr(micara_obj, "disease_col") else "disease"
    unique_vals <- unique(na.omit(micara_obj$metadata[[disease_col]]))

    # Dynamic auto-detection if reference_level is missing
    if (is.null(reference_level)) {
        common_controls <- c("control", "healthy", "ctr", "hc", "wt", "baseline", "normal")
        match_idx <- match(tolower(unique_vals), common_controls)

        if (any(!is.na(match_idx))) {
            reference_level <- unique_vals[which.min(match_idx)]
            if (verbose) {
                message(sprintf("[Info] Auto-selected '%s' as reference level for %s.", reference_level, disease_col))
            }
        } else {
            stop(sprintf(
                "Could not auto-detect baseline level. Please specify 'reference_level'. Choices: %s",
                paste(paste0("'", unique_vals, "'"), collapse = ", ")
            ))
        }
    }

    # Validate supplied reference_level exists
    if (!reference_level %in% unique_vals) {
        stop(sprintf(
            "reference_level '%s' not found in metadata$%s. Available levels: %s",
            reference_level, disease_col, paste(paste0("'", unique_vals, "'"), collapse = ", ")
        ))
    }


    ## ------------------------------------------------------------------
    ## 1. Validate Input & Dependencies
    ## ------------------------------------------------------------------
    # Class & Package Checks: Confirms micara_obj was created by MiCARA functions and verifies that ANCOMBC and phyloseq are installed.

    if (!inherits(micara_obj, "micara_input")) {
        stop(
            "'micara_obj' must be an object of class 'micara_input' (or 'micara_imputed') ",
            "as returned by validate_inputs() or impute_metadata()."
        )
    }

    if (!requireNamespace("ANCOMBC", quietly = TRUE)) {
        stop("Package 'ANCOMBC' is required for run_diff_abundance() but is not installed.")
    }
    if (!requireNamespace("phyloseq", quietly = TRUE)) {
        stop("Package 'phyloseq' is required for run_diff_abundance() but is not installed.")
    }

    ## Parse requested feature types
    # Parses options like "taxa", "pathways", or "both".
    if (any(feature_type %in% c("all", "both"))) {
        feature_type <- c("taxa", "pathways")
    }

    feature_type <- unique(feature_type)
    valid_types <- c("taxa", "pathways")

    invalid <- setdiff(feature_type, valid_types)
    if (length(invalid) > 0) {
        stop(
            "Invalid feature_type specified: ", paste(invalid, collapse = ", "),
            ". Must be 'taxa', 'pathways', or both."
        )
    }

    ## Filter to features actually present in micara_obj, so it can skip empty or missing tables without crashing.
    available_types <- intersect(feature_type, names(micara_obj))
    available_types <- available_types[vapply(available_types, function(ft) {
        !is.null(micara_obj[[ft]]) && ncol(micara_obj[[ft]]) > 0 && nrow(micara_obj[[ft]]) > 0
    }, logical(1))]

    if (length(available_types) == 0) {
        stop(
            "None of the requested feature types (", paste(feature_type, collapse = ", "),
            ") are available or non-empty in 'micara_obj'."
        )
    }

    if (length(available_types) < length(feature_type) && verbose) {
        missing_ft <- setdiff(feature_type, available_types)
        message(
            "[Notice] Skipping requested feature type(s) not present in micara_obj: ",
            paste(missing_ft, collapse = ", ")
        )
    }

    ## ------------------------------------------------------------------
    ## 2. Internal Execution Engine (.run_single)
    ## ------------------------------------------------------------------
    if (verbose) {
        message("[Info] Setting '", reference_level, "' as the reference level for disease comparisons.")
    }
    .run_single <- function(ft) {
        feature_mat <- micara_obj[[ft]]
        metadata <- micara_obj$metadata

        ## Check abundance scale
        current_scale <- micara_obj$abundance_scale$final[[ft]]
        if (!is.null(current_scale) && current_scale != "counts") {
            warning(
                "ANCOM-BC2 log-linear models expect integer count data, but '", ft,
                "' matrix is currently in '", current_scale, "' scale. Results may be affected."
            )
        }

        ## Coerce matrix values to clean integer count format for phyloseq
        feature_mat <- as.matrix(feature_mat)
        if (!is.numeric(feature_mat)) {
            stop("Feature matrix '", ft, "' must contain numeric values.")
        }

        if (!all(feature_mat %% 1 == 0, na.rm = TRUE)) {
            feature_mat <- round(feature_mat)
        }
        storage.mode(feature_mat) <- "integer"

        ## Reference level check & re-leveling around it
        ## Reference level check & re-leveling
        if (!"disease" %in% names(metadata)) {
            stop("'disease' column is missing from metadata.")
        }

        # sanitizes and cleans disease category labels so they don't break R formula syntax or cause string-mismatch errors during model fitting and downstream calculations
        raw_disease <- as.character(metadata$disease)
        clean_disease <- make.names(raw_disease)

        if (!identical(raw_disease, clean_disease)) {
            if (verbose) {
                message("[Notice] Non-standard characters or spaces found in metadata$disease. Sanitizing labels with `make.names()`.")
            }
            metadata$disease <- clean_disease
            # Also update the user's reference_level string to match the sanitized target
            reference_level <- make.names(reference_level)
        }


        disease_levels <- unique(as.character(metadata$disease))

        if (!reference_level %in% disease_levels) {
            stop(
                "reference_level '", reference_level, "' was not found in metadata$disease.\n",
                "Available levels in your dataset: ", paste(paste0("'", disease_levels, "'"), collapse = ", "), "\n",
                "Please specify 'reference_level = \"<your_baseline>\"' in run_diff_abundance()."
            )
        }

        if (verbose) {
            message("[Info] Re-leveling 'disease' with reference level: '", reference_level, "'")
        }
        # set the reference group before building the phyloseq object
        metadata$disease <- stats::relevel(factor(metadata$disease), ref = reference_level)

        ## Determine Fixed Effects
        cur_fix_effects <- fix_effects
        if (is.null(cur_fix_effects)) {
            cur_fix_effects <- setdiff(micara_obj$confounders, c("study_name", "disease"))

            # Automatically prefer continuous 'age' over 'age_category'
            if ("age" %in% cur_fix_effects && "age_category" %in% cur_fix_effects) {
                cur_fix_effects <- setdiff(cur_fix_effects, "age_category")
                if (verbose) {
                    message("[Info] Both 'age' and 'age_category' present. Retaining continuous 'age' and dropping 'age_category' to prevent collinearity.")
                }
            }
        }

        if (!is.null(exclude_covariates)) {
            cur_fix_effects <- setdiff(cur_fix_effects, exclude_covariates)
        }


        ## Validate presence and variance of fixed effect variables and automatically dropping covariates that have zero variance or missing values.
        valid_fix <- character(0)
        for (cov in cur_fix_effects) {
            if (!cov %in% names(metadata)) {
                warning("Covariate '", cov, "' specified in fixed effects was not found in metadata; dropping.")
                next
            }

            col_vals <- metadata[[cov]]
            if (any(is.na(col_vals))) {
                warning(
                    "Covariate '", cov, "' contains missing values. Consider running ",
                    "impute_metadata() first to prevent sample loss during model fitting."
                )
            }

            if (length(unique(stats::na.omit(col_vals))) <= 1) {
                warning("Covariate '", cov, "' has zero variance (constant); dropping from fixed effects.")
                next
            }

            valid_fix <- c(valid_fix, cov)
        }

        fix_formula <- paste(c(valid_fix, "disease"), collapse = " + ")

        ## Determine Random Effects (Multi-Cohort Safeguard)
        # if your dataset only contains a single study, it safely turns off (1 | study_name) to prevent model fitting errors.
        cur_rand_formula <- NULL
        if (!is.null(rand_effect)) {
            if (!rand_effect %in% names(metadata)) {
                if (verbose) {
                    message("[Info] Random effect column '", rand_effect, "' not found in metadata; fitting without random effects.")
                }
            } else {
                n_studies <- length(unique(stats::na.omit(metadata[[rand_effect]])))
                if (n_studies <= 1) {
                    if (verbose) {
                        message(
                            "[Info] '", rand_effect, "' has only ", n_studies,
                            " unique level; fitting fixed-effects model without random intercepts."
                        )
                    }
                } else {
                    cur_rand_formula <- paste0("(1 | ", rand_effect, ")")
                }
            }
        }

        ## Enforce pseudo_sens = TRUE if robust flags are requested
        cur_pseudo_sens <- pseudo_sens
        if (use_robust && !cur_pseudo_sens) {
            if (verbose) {
                message("[Info] Setting 'pseudo_sens = TRUE' to enable robust sensitivity flag calculation.")
            }
            cur_pseudo_sens <- TRUE
        }

        ## ------------------------------------------------------------------
        ## Build Phyloseq Object (Auto-Detect Sample Orientation)
        ## ------------------------------------------------------------------
        rownames(metadata) <- as.character(metadata$sample_id)

        if (all(rownames(metadata) %in% rownames(feature_mat))) {
            # Samples are in ROWS of feature_mat
            feature_mat <- feature_mat[rownames(metadata), , drop = FALSE]
            ps <- phyloseq::phyloseq(
                phyloseq::otu_table(feature_mat, taxa_are_rows = FALSE),
                phyloseq::sample_data(metadata)
            )
        } else if (all(rownames(metadata) %in% colnames(feature_mat))) {
            # Samples are in COLUMNS of feature_mat
            feature_mat <- feature_mat[, rownames(metadata), drop = FALSE]
            ps <- phyloseq::phyloseq(
                phyloseq::otu_table(feature_mat, taxa_are_rows = TRUE),
                phyloseq::sample_data(metadata)
            )
        } else {
            stop(
                "Internal alignment mismatch ", ft, " matrix sample names do not match ",
                "metadata sample IDs. Was micara_obj modified outside MiCARA functions?",
                call. = FALSE
            )
        }

        ## Run ANCOM-BC2
        # executes ANCOMBC::ancombc2() with bias correction, prevalence filtering, and sensitivity analysis.

        if (verbose) {
            rand_form_str <- if (is.null(cur_rand_formula)) "(none)" else cur_rand_formula

            ancom_report <- c(
                paste0("\n========== Running ANCOM-BC2 (", ft, ") =========="),
                paste0("Fixed formula:   ", fix_formula),
                paste0("Random formula:  ", rand_form_str),
                paste0("Reference level: ", reference_level),
                paste0("Method / Alpha:  ", p_adj_method, " / ", alpha),
                "---------------------------------------------------"
            )

            message(paste(ancom_report, collapse = "\n"))
        }

        # Determine current feature type for dynamic message formatting ("taxa" or "pathways")
        cur_type <- if (exists("ft")) ft else "features"

        fit <- withCallingHandlers(
            {
                ANCOMBC::ancombc2(
                    data         = ps,
                    fix_formula  = fix_formula,
                    rand_formula = cur_rand_formula,
                    p_adj_method = p_adj_method,
                    group        = disease_col,
                    prv_cut      = prv_cut,
                    pseudo_sens  = cur_pseudo_sens,
                    struc_zero   = struc_zero,
                    neg_lb       = neg_lb,
                    alpha        = alpha,
                    pairwise     = pairwise,
                    n_cl         = n_cores,
                    verbose      = FALSE
                )
            },
            warning = function(w) {
                # Catch the bias estimation warning from ANCOM-BC2
                if (grepl("The number of taxa used for estimating sample-specific biases", w$message)) {
                    # Replace "taxa" with "pathways" (or current feature type) dynamically
                    custom_msg <- gsub("number of taxa", paste("number of", cur_type), w$message)
                    warning(custom_msg, call. = FALSE)
                    invokeRestart("muffleWarning") # Prevent original warning from printing
                }
            }
        )

        res <- fit$res

        ## Extract Significant Results into Tidy Long Format
        flag_prefix <- if (use_robust) "diff_robust_disease" else "diff_disease"
        diff_cols <- grep(paste0("^", flag_prefix), colnames(res), value = TRUE)

        if (length(diff_cols) == 0 && use_robust) {
            warning("Robust flags ('", flag_prefix, "') not found in ANCOMBC output; falling back to standard differential flags.")
            flag_prefix <- "diff_disease"
            diff_cols <- grep(paste0("^", flag_prefix), colnames(res), value = TRUE)
        }

        comp_levels <- sub(flag_prefix, "", diff_cols)
        feature_col <- if ("taxon" %in% names(res)) "taxon" else names(res)[1]

        significant_list <- lapply(comp_levels, function(dl) {
            diff_col <- paste0(flag_prefix, dl)
            lfc_col <- paste0("lfc_disease", dl)
            se_col <- paste0("se_disease", dl)
            p_col <- paste0("p_disease", dl)
            q_col <- paste0("q_disease", dl)

            if (!diff_col %in% names(res)) {
                return(NULL)
            }

            sig_rows <- res[[diff_col]] & !is.na(res[[diff_col]])

            if (!any(sig_rows)) {
                return(NULL)
            }

            data.frame(
                feature = res[[feature_col]][sig_rows],
                disease = dl,
                lfc = res[[lfc_col]][sig_rows],
                se = res[[se_col]][sig_rows],
                p = res[[p_col]][sig_rows],
                q = res[[q_col]][sig_rows],
                direction = ifelse(res[[lfc_col]][sig_rows] > 0, "enriched", "depleted"),
                stringsAsFactors = FALSE
            )
        })

        significant <- do.call(rbind, significant_list)

        if (is.null(significant) || nrow(significant) == 0) {
            significant <- data.frame(
                feature = character(0),
                disease = character(0),
                lfc = numeric(0),
                se = numeric(0),
                p = numeric(0),
                q = numeric(0),
                direction = character(0),
                stringsAsFactors = FALSE
            )
        } else {
            rownames(significant) <- NULL
        }

        if (verbose) {
            flag_type <- if (use_robust) "robust" else "standard"
            header <- paste0("Significant ", ft, " identified per disease comparison (", flag_type, " flags):")

            if (nrow(significant) > 0) {
                tbl <- table(significant$disease)
                body <- paste(paste0("  ", names(tbl), ": ", tbl), collapse = "\n")
            } else {
                body <- paste0("  None detected at alpha = ", alpha)
            }

            sig_report <- c(
                header,
                body,
                "---------------------------------------------------"
            )

            message(paste(sig_report, collapse = "\n"))
        }

        structure(
            list(
                ancombc2_fit    = fit,
                res             = res,
                significant     = significant,
                feature_type    = ft,
                fix_formula     = fix_formula,
                rand_formula    = cur_rand_formula,
                reference_level = reference_level,
                use_robust      = use_robust
            ),
            class = "micara_diffab"
        )
    }

    ## ------------------------------------------------------------------
    ## 3. Execute Across Feature Types (taxa first then pathways) & Package Output
    ## ------------------------------------------------------------------

    results <- lapply(available_types, .run_single)
    names(results) <- available_types

    if (length(results) == 1) {
        return(results[[1]])
    } else {
        class(results) <- c("micara_diffab_list", "list")
        return(results)
    }
} # If you analyzed one feature type, it returns a single "micara_diffab" object. If you analyzed both, it returns a named list (diff_results$taxa and diff_results$pathways) under class "micara_diffab_list".
