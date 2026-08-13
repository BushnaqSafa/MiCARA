# ================================================================
# MetaPhlAn / HUMAnN Pre-processing Helper (Idempotent)
# ================================================================
clean_shotgun_input <- function(df) {
    if (!is.data.frame(df) || ncol(df) == 0) {
        return(df)
    }

    # Clean Column Headers (MetaPhlAn / HUMAnN Comment Markers)
    first_col_name <- colnames(df)[1]
    if (grepl("^#", first_col_name)) {
        colnames(df)[1] <- gsub("^#\\s*", "", first_col_name)
    }

    # Clean sample column suffixes
    colnames(df) <- gsub("_Abundance$|_Profile$|_pathabundance$", "", colnames(df))

    # Process Feature Column if character or factor
    if (is.character(df[[1]]) || is.factor(df[[1]])) {
        col_vals <- as.character(df[[1]])

        # Filter out unmapped / unintegrated
        unwanted <- grepl("UNINTEGRATED|UNMAPPED", col_vals, ignore.case = TRUE)
        if (any(unwanted)) {
            df <- df[!unwanted, , drop = FALSE]
            col_vals <- as.character(df[[1]])
        }
        # Filter out stratified pathways if raw HUMAnN
        is_stratified <- grepl("\\|", col_vals) &
            grepl("PWY|PATHWAY|UNMAPPED", col_vals, ignore.case = TRUE)

        if (any(is_stratified)) {
            df <- df[!is_stratified, , drop = FALSE]
            col_vals <- as.character(df[[1]])
        }

        # Extract species rank from MetaPhlAn pipe strings
        if (any(grepl("\\|", col_vals))) {
            col_vals <- vapply(strsplit(col_vals, "\\|"), tail, character(1), n = 1L)
        }

        df[[1]] <- col_vals
    }

    return(df)
}

#' Validate and preprocess MiCARA input tables
#'
#' Performs quality control on taxonomic and functional pathway abundance
#' tables and accompanying sample metadata.
#'
#' @param taxa Taxonomic abundance table. Must be a data.frame.
#' @param pathways Functional pathway abundance table. Must be a data.frame.
#' @param metadata Sample metadata. Must contain sample IDs and disease labels.
#' @param disease_col Name of the disease column in metadata. Default "disease".
#' @param disease_min_samples Minimum required samples per disease group. Default 10.
#' @param metadata_missing_cutoff Allowed missing value fraction for covariates. Default 0.30.
#' @param normalise Logical indicating whether feature normalization should be performed.
#' @param verbose Logical; print QC summary. Default TRUE.
#'
#' @return An object of class \code{"micara_input"}.
#' @importFrom utils head tail
#' @export
#' @examples
#' set.seed(123)
#' sample_ids <- paste0("Sample_", 1:20)
#'
#' taxa_mat <- as.data.frame(matrix(
#'     rpois(100, lambda = 50),
#'     nrow = 5, ncol = 20,
#'     dimnames = list(paste0("Taxon_", 1:5), sample_ids)
#' ))
#'
#' path_mat <- as.data.frame(matrix(
#'     rpois(100, lambda = 100),
#'     nrow = 5, ncol = 20,
#'     dimnames = list(paste0("Pathway_", 1:5), sample_ids)
#' ))
#'
#' meta_df <- data.frame(
#'     sample_id = sample_ids,
#'     disease = rep(c("Control", "Disease"), each = 10),
#'     age = rnorm(20, mean = 40, sd = 5),
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
validate_inputs <- function(
  taxa,
  pathways,
  metadata,
  disease_col = "disease",
  disease_min_samples = 10,
  metadata_missing_cutoff = 0.30,
  normalise = TRUE,
  verbose = TRUE
) {
    # Auto-convert matrices, S4 DataFrames, and matrix-like objects to data.frame
    taxa <- tryCatch(as.data.frame(taxa), error = function(e) taxa)
    pathways <- tryCatch(as.data.frame(pathways), error = function(e) pathways)
    metadata <- tryCatch(as.data.frame(metadata), error = function(e) metadata)

    ## ================================================================
    ## 1. Strict Input Type Enforcement
    ## ================================================================
    if (!is.data.frame(taxa)) {
        stop("'taxa' must be a data.frame.", call. = FALSE)
    }
    if (!is.data.frame(pathways)) {
        stop("'pathways' must be a data.frame.", call. = FALSE)
    }
    if (!is.data.frame(metadata)) {
        stop("'metadata' must be a data.frame.", call. = FALSE)
    }

    # Step 0: Preprocess MetaPhlAn/HUMAnN formats non-destructively
    taxa <- clean_shotgun_input(taxa)
    pathways <- clean_shotgun_input(pathways)

    ## ================================================================
    ## 2. Metadata Column & Read Count Aliases Check
    ## ================================================================
    if (!is.character(disease_col) || length(disease_col) != 1L) {
        stop("'disease_col' must be a single column name.", call. = FALSE)
    }

    if (!disease_col %in% names(metadata)) {
        stop("Disease column '", disease_col, "' was not found in metadata.", call. = FALSE)
    }

    if (disease_col != "disease") {
        metadata$disease <- metadata[[disease_col]]
    }

    read_aliases <- c("read_count", "reads", "n_reads", "total_reads", "sequencing_depth", "depth", "number_of_reads")
    if (!"number_reads" %in% names(metadata)) {
        found_alias <- intersect(read_aliases, tolower(names(metadata)))[1]
        if (!is.na(found_alias)) {
            actual_col <- names(metadata)[tolower(names(metadata)) == found_alias][1]
            metadata$number_reads <- metadata[[actual_col]]
        }
    }

    # Auto-populate sample_id from rownames if the column is missing
    if (!"sample_id" %in% colnames(metadata)) {
        if (!is.null(rownames(metadata)) && !all(rownames(metadata) == as.character(seq_len(nrow(metadata))))) {
            metadata$sample_id <- rownames(metadata)
        } else {
            stop("Metadata must contain a 'sample_id' column or valid sample row names.", call. = FALSE)
        }
    }

    if (!"sample_id" %in% names(metadata)) {
        stop("Metadata must contain a column named 'sample_id'.", call. = FALSE)
    }

    ## Universal Metadata Cleaning
    na_strings <- c("missing: not collected", "not collected", "na", "n/a", "", "unknown", "null", "nan", "none")
    metadata[] <- lapply(metadata, function(col) {
        if (is.character(col) || is.factor(col)) {
            col_char <- trimws(as.character(col))
            col_clean <- tolower(col_char)
            col_char[col_clean %in% na_strings] <- NA_character_
            col <- col_char
        }
        if (is.character(col)) {
            num_converted <- suppressWarnings(as.numeric(col))
            if (sum(!is.na(col)) > 0 && sum(!is.na(num_converted)) == sum(!is.na(col))) {
                return(num_converted)
            }
        }
        return(col)
    })

    ## ================================================================
    ## 3. Convert Abundance Data Frames & Separate Feature Name Column
    ## ================================================================
    prepare_abundance_table <- function(df, label) {
        if (ncol(df) < 1L) {
            stop("'", label, "' must contain at least one column.", call. = FALSE)
        }

        col1_name <- colnames(df)[1]
        col1_is_feature <- is.character(df[[1]]) || is.factor(df[[1]]) ||
            tolower(col1_name) %in% c("feature_id", "feature", "taxon", "taxa", "pathway", "gene", "clade", "id", "#otu id", "otu_id") ||
            grepl("^#", col1_name)

        if (col1_is_feature) {
            if (ncol(df) < 2L) {
                stop("'", label, "' must contain a feature-name column and at least one sample column.", call. = FALSE)
            }
            feature_names <- as.character(df[[1]])
            mat_df <- df[, -1, drop = FALSE]
        } else {
            feature_names <- rownames(df)
            if (is.null(feature_names) || identical(feature_names, as.character(seq_len(nrow(df))))) {
                feature_names <- paste0("feature_", seq_len(nrow(df)))
            }
            mat_df <- df
        }

        if (anyNA(feature_names) || any(feature_names == "")) {
            stop("'", label, "' contains missing or empty feature names.", call. = FALSE)
        }

        if (anyDuplicated(feature_names)) {
            dup_names <- unique(feature_names[duplicated(feature_names)])
            stop("'", label, "' contains duplicated feature names: ", paste(head(dup_names, 5L), collapse = ", "), call. = FALSE)
        }

        non_num <- !vapply(mat_df, is.numeric, logical(1))
        if (any(non_num)) {
            stop("'", label, "' contains non-numeric columns.", call. = FALSE)
        }

        mat <- as.matrix(mat_df)
        storage.mode(mat) <- "double"
        rownames(mat) <- feature_names

        if (is.null(colnames(mat)) || anyNA(colnames(mat)) || any(colnames(mat) == "")) {
            stop("'", label, "' contains missing or empty sample IDs.", call. = FALSE)
        }

        if (anyDuplicated(colnames(mat))) {
            stop("'", label, "' contains duplicated sample IDs.", call. = FALSE)
        }

        if (anyNA(mat) || any(!is.finite(mat))) {
            stop("'", label, "' contains missing or non-finite abundance values.", call. = FALSE)
        }

        if (any(mat < 0)) {
            stop("'", label, "' contains negative abundance values.", call. = FALSE)
        }

        return(mat)
    }

    taxa <- prepare_abundance_table(taxa, "taxa")
    pathways <- prepare_abundance_table(pathways, "pathways")

    ## ================================================================
    ## 4. Sample ID Checks & Alignment
    ## ================================================================
    taxa_samples <- colnames(taxa)
    pathway_samples <- colnames(pathways)
    metadata_samples <- as.character(metadata$sample_id)

    if (anyNA(metadata_samples) || any(metadata_samples == "")) {
        stop("'metadata$sample_id' contains missing or empty sample IDs.", call. = FALSE)
    }

    if (anyDuplicated(metadata_samples)) {
        stop("'metadata$sample_id' contains duplicated sample IDs.", call. = FALSE)
    }

    if (!identical(taxa_samples, pathway_samples)) {
        stop("Taxa and pathway sample IDs are not identical and in the same order.", call. = FALSE)
    }

    missing_metadata <- setdiff(taxa_samples, metadata_samples)
    if (length(missing_metadata) > 0L) {
        stop(
            length(missing_metadata),
            " sample(s) present in taxa/pathways but absent from metadata: ",
            paste(head(missing_metadata, 5L), collapse = ", "),
            call. = FALSE
        )
    }

    extra_metadata <- setdiff(metadata_samples, taxa_samples)
    if (length(extra_metadata) > 0L && verbose) {
        warning(length(extra_metadata), " metadata sample(s) are not present in abundance tables and will be removed.", call. = FALSE)
    }

    ## Align all datasets
    common_samples <- taxa_samples
    taxa <- taxa[, common_samples, drop = FALSE]
    pathways <- pathways[, common_samples, drop = FALSE]
    metadata <- metadata[match(common_samples, metadata$sample_id), , drop = FALSE]

    ## ================================================================
    ## 5. Disease Group Filtering
    ## ================================================================
    if (anyNA(metadata$disease) || any(metadata$disease == "")) {
        warning("Some samples have missing disease labels and will be removed.", call. = FALSE)
        keep <- !is.na(metadata$disease) & metadata$disease != ""
        metadata <- metadata[keep, , drop = FALSE]
        taxa <- taxa[, keep, drop = FALSE]
        pathways <- pathways[, keep, drop = FALSE]
    }

    disease_counts <- table(metadata$disease)
    keep_diseases <- names(disease_counts[disease_counts >= disease_min_samples])
    removed_diseases <- setdiff(names(disease_counts), keep_diseases)

    keep_sample_ids <- metadata$sample_id[metadata$disease %in% keep_diseases]
    taxa <- taxa[, keep_sample_ids, drop = FALSE]
    pathways <- pathways[, keep_sample_ids, drop = FALSE]
    metadata <- metadata[match(keep_sample_ids, metadata$sample_id), , drop = FALSE]

    ## ================================================================
    ## 6. Zero-Sum Filtering
    ## ================================================================
    taxa_zero_rows <- rowSums(taxa) == 0
    if (any(taxa_zero_rows)) {
        if (verbose) warning(sprintf("Removing %d all-zero taxa feature(s).", sum(taxa_zero_rows)), call. = FALSE)
        taxa <- taxa[!taxa_zero_rows, , drop = FALSE]
    }

    path_zero_rows <- rowSums(pathways) == 0
    if (any(path_zero_rows)) {
        if (verbose) warning(sprintf("Removing %d all-zero pathway feature(s).", sum(path_zero_rows)), call. = FALSE)
        pathways <- pathways[!path_zero_rows, , drop = FALSE]
    }

    valid_samples <- (colSums(taxa) > 0) & (colSums(pathways) > 0)
    if (any(!valid_samples)) {
        dropped_ids <- colnames(taxa)[!valid_samples]
        if (verbose) warning(sprintf("Removing %d sample(s) with depth < 1 reads: %s", sum(!valid_samples), paste(dropped_ids, collapse = ", ")), call. = FALSE)
        taxa <- taxa[, valid_samples, drop = FALSE]
        pathways <- pathways[, valid_samples, drop = FALSE]
        metadata <- metadata[match(colnames(taxa), metadata$sample_id), , drop = FALSE]
    }

    ## ================================================================
    ## 7. Case-Insensitive Covariate / Confounder Assessment
    ## ================================================================
    recognised_covariates <- c("study_name", "age", "age_category", "gender", "sex", "bmi", "body_site", "site")
    meta_cols <- names(metadata)
    is_covariate <- tolower(meta_cols) %in% tolower(recognised_covariates) &
        !tolower(meta_cols) %in% c("sample_id", "disease", tolower(disease_col), "number_reads")
    available_covariates <- meta_cols[is_covariate]

    if (length(available_covariates) > 0L) {
        missing_fraction <- vapply(available_covariates, function(x) mean(is.na(metadata[[x]])), numeric(1))
        retained_covariates <- available_covariates[missing_fraction <= metadata_missing_cutoff]
        removed_covariates <- available_covariates[missing_fraction > metadata_missing_cutoff]
    } else {
        retained_covariates <- character(0)
        removed_covariates <- character(0)
    }

    ## ================================================================
    ## 8. Abundance Scale Detection & Re-construction
    ## ================================================================
    detect_scale <- function(mat) {
        if (is.data.frame(mat)) {
            num_cols <- vapply(mat, is.numeric, logical(1))
            mat <- as.matrix(mat[, num_cols, drop = FALSE])
        } else {
            mat <- as.matrix(mat)
        }

        if (length(mat) == 0) {
            stop("Unrecognised abundance scale. Unable to determine scale.", call. = FALSE)
        }

        # 1. Negative values -> log/transformed
        if (any(mat < 0, na.rm = TRUE)) {
            return("transformed/log")
        }

        # 2. Integer counts
        is_integer_like <- all(abs(mat - round(mat)) < 1e-5, na.rm = TRUE)
        if (is_integer_like) {
            return("counts")
        }

        totals <- colSums(mat, na.rm = TRUE)
        median_total <- stats::median(totals, na.rm = TRUE)
        max_val <- max(mat, na.rm = TRUE)
        max_total <- max(totals, na.rm = TRUE)

        # 3. Proportion scale (0 to 1)
        # Max value <= 1.05 & max total <= 1.05 (handles filtered pathways smoothly)
        if (max_val <= 1.05 && max_total <= 1.05) {
            return("proportion")
        }

        # 4. Percentage scale (0 to 100)
        # Max value <= 100.05 & sample sums average ~100
        if (max_val <= 100.05 && abs(median_total - 100) <= 10) {
            return("percentage")
        }

        # 5. Fallback: Throw informative error for ambiguous continuous scale
        stop("Unrecognised abundance scale. Unable to determine scale.", call. = FALSE)
    }

    taxa_scale <- detect_scale(taxa)
    pathways_scale <- detect_scale(pathways)

    has_reads <- "number_reads" %in% names(metadata) && !all(is.na(metadata$number_reads)) && is.numeric(metadata$number_reads)

    reconstruct_counts <- function(mat, scale, sample_reads) {
        if (scale == "counts") {
            return(mat)
        }
        prop_mat <- if (scale == "percentage") mat / 100 else mat
        round(sweep(prop_mat, 2, sample_reads, FUN = "*"))
    }

    process_dataset <- function(mat, scale, name) {
        if (scale %in% c("percentage", "proportion")) {
            if (has_reads) {
                message("[", name, "] Back-calculating estimated raw counts using 'metadata$number_reads'.")
                return(reconstruct_counts(mat, scale, metadata$number_reads))
            } else {
                warning("[", name, "] Input is in '", scale, "' scale and 'metadata$number_reads' was not provided.", call. = FALSE)
            }
        }
        return(mat)
    }

    taxa <- process_dataset(taxa, taxa_scale, "taxa")
    pathways <- process_dataset(pathways, pathways_scale, "pathways")

    final_taxa_scale <- if (has_reads && taxa_scale %in% c("percentage", "proportion")) "counts" else taxa_scale
    final_pathways_scale <- if (has_reads && pathways_scale %in% c("percentage", "proportion")) "counts" else pathways_scale

    ## ================================================================
    ## 9. QC Report Output
    ## ================================================================
    if (verbose) {
        removed_dis_str <- if (length(removed_diseases) > 0L) paste(removed_diseases, collapse = ", ") else "None"
        retained_cov_str <- if (length(retained_covariates) > 0L) paste(retained_covariates, collapse = ", ") else "None"
        removed_cov_str <- if (length(removed_covariates) > 0L) paste(removed_covariates, collapse = ", ") else "None"

        message(paste(c(
            "\n========== MiCARA Input QC ==========",
            paste0("Samples retained:        ", ncol(taxa)),
            paste0("Diseases retained:       ", length(keep_diseases)),
            paste0("Diseases removed:        ", removed_dis_str),
            paste0("Taxa scale detected:     ", taxa_scale),
            paste0("Pathway scale detected:  ", pathways_scale),
            paste0("Covariates retained:     ", retained_cov_str),
            paste0("Covariates removed:      ", removed_cov_str),
            "====================================="
        ), collapse = "\n"))
    }

    ## ================================================================
    ## 10. Return Structured Object
    ## ================================================================
    structure(
        list(
            taxa = as.data.frame(taxa, check.names = FALSE),
            pathways = as.data.frame(pathways, check.names = FALSE),
            metadata = metadata,
            abundance_type = list(
                detected = list(
                    taxa = taxa_scale,
                    pathways = pathways_scale
                ),
                final = list(
                    taxa = final_taxa_scale,
                    pathways = final_pathways_scale
                )
            ),
            confounders = retained_covariates,
            removed_covariates = removed_covariates,
            removed_diseases = removed_diseases
        ),
        class = "micara_input"
    )
}
