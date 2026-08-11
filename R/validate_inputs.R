
#' Validate and preprocess MiCARA input tables
#'
#' Performs quality control on taxonomic and functional pathway abundance
#' tables and accompanying sample metadata.
#'
#' @param taxa Taxonomic abundance table. Features must be rows and samples
#'   columns. The first column may contain feature names.
#' @param pathways Functional pathway abundance table. Features must be rows
#'   and samples columns. The first column may contain feature names.
#' @param metadata Sample metadata. Must contain sample IDs and disease labels.
#' @param disease_col Name of the disease column in metadata.
#' @param disease_min_samples Minimum number of samples required to retain a
#'   disease group. Default is 10.
#' @param metadata_missing_cutoff Maximum allowed fraction of missing values
#'   for an optional covariate. Default is 0.30.
#' @param normalise Logical indicating whether feature normalization should be performed. Default is TRUE.
#' @param verbose Logical; print a QC summary. Default TRUE.
#'
#' @return An object of class \code{"micara_input"}.
#'
#' @export
validate_inputs <- function(
    taxa,
    pathways,
    metadata,
    disease_col             = "disease",
    disease_min_samples     = 10,
    metadata_missing_cutoff = 0.30,
    normalise                = TRUE,
    verbose                  = TRUE
) {
  
  ## ================================================================
  ## 1. Basic input validation (Object type and auto-coercion)
  ## ================================================================
  
  coerce_to_df <- function(x, name) {
    
    if (is.matrix(x) ||
        is.data.frame(x) ||
        inherits(x, c("DFrame", "DataFrame", "data.table"))) {
      
      return(as.data.frame(x, stringsAsFactors = FALSE))
    }
    
    stop(
      sprintf(
        "'%s' must be a matrix, data.frame, tibble, or Bioconductor DataFrame (received class: %s)",
        name,
        paste(class(x), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  taxa     <- coerce_to_df(taxa, "taxa")
  pathways <- coerce_to_df(pathways, "pathways")
  metadata <- coerce_to_df(metadata, "metadata")
  
  if (!is.character(disease_col) || length(disease_col) != 1L) {
    stop("'disease_col' must be a single column name.", call. = FALSE)
  }
  
  if (!disease_col %in% names(metadata)) {
    stop(
      "Disease column '", disease_col,
      "' was not found in metadata.",
      call. = FALSE
    )
  }
  
  if (!"sample_id" %in% names(metadata)) {
    stop(
      "Metadata must contain a column named 'sample_id'.",
      call. = FALSE
    )
  }
  
  ## ================================================================
  ## Universal Metadata Cleaning (Convert common missing-value strings to NA and auto-detect
  ## numeric columns)
  ## ================================================================
  na_strings <- c(
    "missing: not collected",
    "not collected",
    "na",
    "n/a",
    "",
    "unknown",
    "null",
    "nan",
    "none"
  )
  
  metadata[] <- lapply(metadata, function(col) {
    # 1. Convert string representation of missing values to real R NAs
    if (is.character(col) || is.factor(col)) {
      
      col_char <- trimws(as.character(col))
      
      ## Case-insensitive matching
      col_clean <- tolower(col_char)
      
      col_char[col_clean %in% na_strings] <- NA_character_
      
      col <- col_char
    }
    
    # 2. Auto-detect if a character column is actually numeric
    if (is.character(col)) {
      
      num_converted <- suppressWarnings(as.numeric(col))
      
      ## Convert to numeric only when every non-missing value
      ## can be converted successfully
      if (
        sum(!is.na(col)) > 0 &&
        sum(!is.na(num_converted)) == sum(!is.na(col))
      ) {
        return(num_converted)
      }
    }
    
    return(col)
  })
  
  ## ================================================================
  ## 2. Convert feature tables to feature x sample matrices
  ##
  ## CAMDA format:
  ##
  ##   feature_name | sample1 | sample2 | ...
  ##
  ## After this step:
  ##
  ##   rownames = feature names
  ##   columns  = sample IDs
  ## ================================================================
  prepare_abundance_table <- function(df, label) {
    
    if (ncol(df) < 1L) {
      stop(
        "'", label, "' must contain at least one sample column.",
        call. = FALSE
      )
    }
    
    ## --------------------------------------------------------------
    ## Determine whether feature names are stored in rownames
    ## or in the first column.
    ## --------------------------------------------------------------
    
    all_numeric <- all(vapply(df, is.numeric, logical(1)))
    
    has_custom_rownames <-
      !is.null(rownames(df)) &&
      length(rownames(df)) == nrow(df) &&
      all(!is.na(rownames(df))) &&
      all(rownames(df) != "") &&
      !identical(rownames(df), as.character(seq_len(nrow(df))))
    
    ## ==============================================================
    ## FORMAT A: Feature names are already stored in rownames
    ## ==============================================================
    
    if (all_numeric && has_custom_rownames) {
      
      feature_names <- rownames(df)
      mat_df <- df
      
    } else {
      
      ## ============================================================
      ## FORMAT B: First column contains feature names
      ## ============================================================
      
      if (ncol(df) < 2L) {
        stop(
          "'", label,
          "' must contain a feature-name column and at least one sample column.",
          call. = FALSE
        )
      }
      
      feature_names <- as.character(df[[1]])
      mat_df <- df[, -1, drop = FALSE]
      
      ## All abundance columns must be numeric
      numeric_remaining <- vapply(
        mat_df,
        is.numeric,
        logical(1)
      )
      
      if (!all(numeric_remaining)) {
        stop(
          "'", label,
          "' contains non-numeric sample columns. ",
          "Expected a feature-name column followed by numeric abundance values.",
          call. = FALSE
        )
      }
    }
    
    ## --------------------------------------------------------------
    ## Validate feature names
    ## --------------------------------------------------------------
    
    if (anyNA(feature_names) || any(feature_names == "")) {
      stop(
        "'", label, "' contains missing or empty feature names.",
        call. = FALSE
      )
    }
    
    if (anyDuplicated(feature_names)) {
      
      duplicated_names <- unique(
        feature_names[duplicated(feature_names)]
      )
      
      stop(
        "'", label,
        "' contains duplicated feature names: ",
        paste(head(duplicated_names, 5L), collapse = ", "),
        if (length(duplicated_names) > 5L) " ..." else "",
        call. = FALSE
      )
    }
    
    ## --------------------------------------------------------------
    ## Convert abundance table to numeric matrix
    ## --------------------------------------------------------------
    
    mat <- as.matrix(mat_df)
    storage.mode(mat) <- "double"
    
    ## Assign feature names
    rownames(mat) <- feature_names
    
    ## --------------------------------------------------------------
    ## Validate sample IDs
    ## --------------------------------------------------------------
    
    if (
      is.null(colnames(mat)) ||
      anyNA(colnames(mat)) ||
      any(colnames(mat) == "")
    ) {
      
      stop(
        "'", label, "' contains missing or empty sample IDs.",
        call. = FALSE
      )
    }
    
    if (anyDuplicated(colnames(mat))) {
      
      duplicated_ids <- unique(
        colnames(mat)[duplicated(colnames(mat))]
      )
      
      stop(
        "'", label,
        "' contains duplicated sample IDs: ",
        paste(head(duplicated_ids, 5L), collapse = ", "),
        if (length(duplicated_ids) > 5L) " ..." else "",
        call. = FALSE
      )
    }
    
    ## --------------------------------------------------------------
    ## Validate abundance values
    ## --------------------------------------------------------------
    
    if (anyNA(mat) || any(!is.finite(mat))) {
      stop(
        "'", label,
        "' contains missing or non-finite abundance values.",
        call. = FALSE
      )
    }
    
    if (any(mat < 0)) {
      stop(
        "'", label,
        "' contains negative abundance values.",
        call. = FALSE
      )
    }
    
    mat
  }
  
  ## ================================================================
  ## 3. Sample ID checks
  ## ================================================================
  
  taxa_samples <- colnames(taxa)
  pathway_samples <- colnames(pathways)
  metadata_samples <- as.character(metadata$sample_id)
  
  if (anyNA(metadata_samples) || any(metadata_samples == "")) {
    stop(
      "'metadata$sample_id' contains missing or empty sample IDs.",
      call. = FALSE
    )
  }
  
  if (anyDuplicated(metadata_samples)) {
    duplicated_ids <- unique(
      metadata_samples[duplicated(metadata_samples)]
    )
    
    stop(
      "'metadata$sample_id' contains duplicated sample IDs: ",
      paste(head(duplicated_ids, 5L), collapse = ", "),
      if (length(duplicated_ids) > 5L) " ..." else "",
      call. = FALSE
    )
  }
  
  ## Taxa and pathways must contain exactly the same samples
  if (!setequal(taxa_samples, pathway_samples)) {
    
    taxa_only <- setdiff(taxa_samples, pathway_samples)
    pathway_only <- setdiff(pathway_samples, taxa_samples)
    
    stop(
      "Taxa and pathway sample IDs do not match.\n",
      "Samples only in taxa: ",
      paste(head(taxa_only, 5L), collapse = ", "),
      if (length(taxa_only) > 5L) " ..." else "",
      "\nSamples only in pathways: ",
      paste(head(pathway_only, 5L), collapse = ", "),
      if (length(pathway_only) > 5L) " ..." else "",
      call. = FALSE
    )
  }
  
  ## Metadata must contain all abundance samples
  missing_metadata <- setdiff(taxa_samples, metadata_samples)
  
  if (length(missing_metadata) > 0L) {
    
    stop(
      length(missing_metadata),
      " sample(s) present in taxa/pathways but absent from metadata: ",
      paste(head(missing_metadata, 5L), collapse = ", "),
      if (length(missing_metadata) > 5L) " ..." else "",
      call. = FALSE
    )
  }
  
  ## Remove metadata-only samples if present
  extra_metadata <- setdiff(metadata_samples, taxa_samples)
  
  if (length(extra_metadata) > 0L && verbose) {
    
    warning(
      length(extra_metadata),
      " metadata sample(s) are not present in the abundance tables ",
      "and will be removed.",
      call. = FALSE
    )
  }
  
  ## ================================================================
  ## 4. Align all three datasets by sample ID
  ## ================================================================
  
  common_samples <- taxa_samples
  
  taxa <- taxa[, common_samples, drop = FALSE]
  pathways <- pathways[, common_samples, drop = FALSE]
  
  metadata <- metadata[
    match(common_samples, metadata$sample_id),
    ,
    drop = FALSE
  ]
  
  ## Rename disease column internally to 'disease'
  ## This allows users to supply another column name.
  metadata$disease <- metadata[[disease_col]]
  
  ## Final alignment check
  if (!identical(colnames(taxa), colnames(pathways)) ||
      !identical(colnames(taxa), metadata$sample_id)) {
    
    stop(
      "Internal alignment error: taxa, pathways and metadata ",
      "sample IDs are not identical after alignment.",
      call. = FALSE
    )
  }
  
  ## ================================================================
  ## 5. Disease filtering
  ## ================================================================
  
  if (anyNA(metadata$disease) || any(metadata$disease == "")) {
    warning(
      "Some samples have missing disease labels and will be removed.",
      call. = FALSE
    )
    
    keep <- !is.na(metadata$disease) &
      metadata$disease != ""
    
    metadata <- metadata[keep, , drop = FALSE]
    taxa <- taxa[, keep, drop = FALSE]
    pathways <- pathways[, keep, drop = FALSE]
  }
  
  disease_counts <- table(metadata$disease)
  
  keep_diseases <- names(
    disease_counts[disease_counts >= disease_min_samples]
  )
  
  removed_diseases <- setdiff(
    names(disease_counts),
    keep_diseases
  )
  
  ## IMPORTANT:
  ## Select sample IDs rather than applying the metadata logical
  ## vector directly to abundance-table columns.
  keep_sample_ids <- metadata$sample_id[
    metadata$disease %in% keep_diseases
  ]
  
  taxa <- taxa[
    ,
    keep_sample_ids,
    drop = FALSE
  ]
  
  pathways <- pathways[
    ,
    keep_sample_ids,
    drop = FALSE
  ]
  
  metadata <- metadata[
    match(keep_sample_ids, metadata$sample_id),
    ,
    drop = FALSE
  ]
  
  ## Final disease-filter alignment check
  if (!identical(colnames(taxa), colnames(pathways)) ||
      !identical(colnames(taxa), metadata$sample_id)) {
    
    stop(
      "Internal alignment error after disease filtering.",
      call. = FALSE
    )
  }
  
  ## ================================================================
  ## 6. Remove all-zero taxa or pathways and empty samples
  ## ================================================================
  # Step A: Drop features (rows) with 0 total counts across all samples
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
  
  # Step B: Drop samples (columns) with 0 total sequencing depth
  
  
  valid_samples <- (colSums(taxa) > 0) & (colSums(pathways) > 0)
  
  if (any(!valid_samples)) {
    dropped_ids <- colnames(taxa)[!valid_samples]
    if (verbose) {
      warning(sprintf("Removing %d sample(s) with sequencing depth < %d reads: %s", 
                      sum(!valid_samples), 
                      paste(dropped_ids, collapse = ", ")), call. = FALSE)
    }
    
    # Step C: Sync matrices and metadata
    taxa     <- taxa[, valid_samples, drop = FALSE]
    pathways <- pathways[, valid_samples, drop = FALSE]
    metadata <- metadata[match(colnames(taxa), metadata$sample_id), , drop = FALSE]
  }

  ## ================================================================
  ## 7. Optional metadata / confounder assessment
  ## ================================================================
  
  recognised_covariates <- c(
    "study_name",
    "age",
    "age_category",
    "gender",
    "BMI"
  )
  
  available_covariates <- intersect(
    recognised_covariates,
    names(metadata)
  )
  
  ## Warn about expected columns with different capitalisation
  lower_names <- tolower(names(metadata))
  
  possible_mismatches <- recognised_covariates[
    !(recognised_covariates %in% names(metadata)) &
      tolower(recognised_covariates) %in% lower_names
  ]
  
  if (length(possible_mismatches) > 0L) {
    
    actual_names <- names(metadata)[
      match(tolower(possible_mismatches), lower_names)
    ]
    
    warning(
      "Possible metadata naming mismatch: ",
      paste(
        paste0(
          possible_mismatches,
          " found as '",
          actual_names,
          "'"
        ),
        collapse = "; "
      ),
      ". Rename columns to the expected names if they should be used.",
      call. = FALSE
    )
  }
  
  ## Missingness assessment
  missing_fraction <- vapply(
    available_covariates,
    function(x) {
      mean(is.na(metadata[[x]]))
    },
    numeric(1)
  )
  
  retained_covariates <- names(
    missing_fraction[
      missing_fraction <= metadata_missing_cutoff
    ]
  )
  
  removed_covariates <- names(
    missing_fraction[
      missing_fraction > metadata_missing_cutoff
    ]
  )
  
  ## Clean continuous covariates (convert string NAs to real NAs and coerce to numeric)
  for (num_col in c("age", "BMI")) {
    if (num_col %in% names(metadata)) {
      vals <- as.character(metadata[[num_col]])
      vals[vals %in% c("NA", "N/A", "", "Unknown", "unknown", "null")] <- NA
      metadata[[num_col]] <- suppressWarnings(as.numeric(vals))
    }
  }
  ## ------------------------------------------------------------------
  ## 8. HELPER: Detect abundance scale
  ## ------------------------------------------------------------------
  
  detect_scale <- function(mat) {
    if (any(mat < 0, na.rm = TRUE)) {
      return("transformed/log")
    }
    
    totals       <- colSums(mat, na.rm = TRUE)
    median_total <- stats::median(totals, na.rm = TRUE)
    max_val      <- max(mat, na.rm = TRUE)
    
    # Proportion: bounded by 1.0 (holds true even after removing unmapped rows)
    if (max_val <= 1.0 && median_total <= 1.0) {
      return("proportion")
    } 
    
    # Percentage: values <= 100 with sample totals centered near or below 100
    if (max_val <= 100 && (abs(median_total - 100) < 15 || max_val > 1)) {
      return("percentage")
    } 
    
    # Raw counts: large integer abundances
    if (median_total > 1) {
      return("counts")
    }
    
    stop("Unable to determine abundance scale. Median sample total = ", round(median_total, 2))
  }
  
  # Detect scale for both matrices
  taxa_scale     <- detect_scale(taxa)
  pathways_scale <- detect_scale(pathways)
  
  
  ## ------------------------------------------------------------------
  ## 9. HELPER: Back-calculate raw counts from relative abundance 
  ## ------------------------------------------------------------------
  
  reconstruct_counts <- function(mat, scale, sample_reads) {
    if (scale == "counts") {
      return(mat)
    }
    
    if (is.null(sample_reads) || any(is.na(sample_reads))) {
      warning("Cannot reconstruct counts: 'number_reads' column is missing or contains NAs in metadata. Proceeding without reconstruction.")
      return(mat)
    }
    
    if (scale == "percentage") {
      prop_mat <- mat / 100
    } else if (scale == "proportion") {
      prop_mat <- mat
    } else {
      warning("Scale '", scale, "' cannot be automatically converted to counts. Returning original matrix.")
      return(mat)
    }
    # Reconstruct counts: proportion * total sample reads, rounded to nearest integer
    count_mat <- round(sweep(prop_mat, 2, sample_reads, FUN = "*"))
    
    message(
      "\n[CAVEAT] Back-calculated estimated counts by multiplying relative abundances ",
      "by metadata$number_reads. Note that these are approximations, and true raw integer ",
      "counts are always preferred for ANCOM-BC2 models.\n"
    )
    
    return(count_mat)
  } 
  
  ## ------------------------------------------------------------------
  ## 10. VERIFICATION & RECONSTRUCTION PIPELINE
  ## ------------------------------------------------------------------
  # Check if total read counts are available and complete
  has_reads <- "number_reads" %in% names(metadata) && !any(is.na(metadata$number_reads))
  
  process_dataset <- function(mat, scale, name) {
    if (scale %in% c("percentage", "proportion")) {
      if (has_reads) {
        warning(
          "[", name, "] Input is in '", scale, "' scale. ANCOM-BC2 strictly expects raw counts. ",
          "Back-calculating estimated raw counts using 'metadata$number_reads'.",
          call. = FALSE
        )
        return(reconstruct_counts(mat, scale, metadata$number_reads))
      } else {
        warning(
          "[", name, "] Input is in '", scale, "' scale and 'metadata$number_reads' was not provided. ",
          "ANCOM-BC2 strictly expects raw counts, but proceeding with relative abundances as requested.",
          call. = FALSE
        )
      }
    }
    return(mat)
  }
  
  taxa     <- process_dataset(taxa, taxa_scale, "taxa")
  pathways <- process_dataset(pathways, pathways_scale, "pathways")
  

  ## ================================================================
  ## 11. QC report
  ## ================================================================
  
  if (verbose) {
    
    cat("\n")
    cat("========== MiCARA Input QC ==========\n")
    
    cat(
      "Samples retained:       ",
      ncol(taxa),
      "\n"
    )
    
    cat(
      "Diseases retained:      ",
      length(keep_diseases),
      "\n"
    )
    
    cat(
      "Disease minimum n:      ",
      disease_min_samples,
      "\n"
    )
    
    if (length(removed_diseases) > 0L) {
      
      cat(
        "Diseases removed:       ",
        paste(removed_diseases, collapse = ", "),
        "\n"
      )
      
    } else {
      
      cat(
        "Diseases removed:       None\n"
      )
    }
    
    cat(
      "Taxa scale detected:    ",
      taxa_scale,
      "\n"
    )
    
    cat(
      "Pathway scale detected: ",
      pathways_scale,
      "\n"
    )
    
    cat(
      "Normalised:             ",
      normalise,
      "\n"
    )
    
    if (length(retained_covariates) > 0L) {
      
      cat(
        "Covariates retained:    ",
        paste(retained_covariates, collapse = ", "),
        "\n"
      )
      
    } else {
      
      cat(
        "Covariates retained:    None\n"
      )
    }
    
    if (length(removed_covariates) > 0L) {
      
      cat(
        "Covariates removed:     ",
        paste(removed_covariates, collapse = ", "),
        " (> ",
        metadata_missing_cutoff * 100,
        "% missing)\n",
        sep = ""
      )
      
    } else {
      
      cat(
        "Covariates removed:     None\n"
      )
    }
    
    cat(
      "=====================================\n\n"
    )
  }
  
  ## ================================================================
  ## 12. Return MiCARA input object
  ## ================================================================
  
  structure(
    list(
      taxa = taxa,
      pathways = pathways,
      metadata = metadata,
      
      abundance_type = list(
        taxa = taxa_scale,
        pathways = pathways_scale
      ),
      
      confounders = retained_covariates,
      
      removed_covariates = removed_covariates,
      
      removed_diseases = removed_diseases
    ),
    class = "micara_input"
  )
}
