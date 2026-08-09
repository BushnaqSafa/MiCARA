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
#' @param normalise Logical; if TRUE, abundance tables are converted to
#'   percentage relative abundance. Default TRUE.
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
  ## 1. Basic input validation
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
    
    if (ncol(df) < 2L) {
      stop(
        "'", label,
        "' must contain a feature-name column and at least one sample column.",
        call. = FALSE
      )
    }
    
    ## --------------------------------------------------------------
    ## First column is treated as the feature-name column when all
    ## remaining columns are numeric.
    ## --------------------------------------------------------------
    
    numeric_remaining <- vapply(
      df[, -1, drop = FALSE],
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
    
    feature_names <- as.character(df[[1]])
    
    if (anyNA(feature_names) || any(feature_names == "")) {
      stop(
        "'", label,
        "' contains missing or empty feature names.",
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
    
    ## Remove feature-name column
    mat <- df[, -1, drop = FALSE]
    
    ## Convert to numeric matrix
    mat <- as.matrix(mat)
    storage.mode(mat) <- "double"
    
    ## Assign feature names
    rownames(mat) <- feature_names
    
    ## Sample IDs must exist
    if (is.null(colnames(mat)) ||
        anyNA(colnames(mat)) ||
        any(colnames(mat) == "")) {
      
      stop(
        "'", label,
        "' contains missing or empty sample IDs.",
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
    
    ## Check values
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
  
  taxa <- prepare_abundance_table(taxa, "taxa")
  pathways <- prepare_abundance_table(pathways, "pathways")
  
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
  ## 6. Optional metadata / confounder assessment
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
  
  ## ================================================================
  ## 7. Abundance scale detection
  ## ================================================================
  
  detect_scale <- function(mat, label) {
    
    totals <- colSums(mat)
    
    if (any(!is.finite(totals)) || any(totals <= 0)) {
      stop(
        "'", label,
        "' contains sample(s) with zero or invalid total abundance.",
        call. = FALSE
      )
    }
    
    median_total <- median(totals)
    
    if (abs(median_total - 100) <= 5) {
      return("percentage")
    }
    
    if (abs(median_total - 1) <= 0.05) {
      return("proportion")
    }
    
    ## Counts are better identified by integer-like values
    ## and substantially larger sample totals.
    if (
      median_total > 1000 &&
      all(abs(mat - round(mat)) < 1e-8)
    ) {
      return("counts")
    }
    
    stop(
      "Unable to determine abundance scale for '", label,
      "'. Median sample total = ",
      round(median_total, 4),
      ". Expected percentage (~100), proportion (~1), ",
      "or count data (>1000 with integer-like values).",
      call. = FALSE
    )
  }
  
  taxa_scale <- detect_scale(taxa, "taxa")
  pathways_scale <- detect_scale(pathways, "pathways")
  
  ## ================================================================
  ## 8. Normalise to percentage scale
  ## ================================================================
  
  normalize_to_percentage <- function(mat, scale) {
    
    if (scale == "percentage") {
      return(mat)
    }
    
    if (scale == "proportion") {
      return(mat * 100)
    }
    
    if (scale == "counts") {
      
      totals <- colSums(mat)
      
      return(
        sweep(
          mat,
          2,
          totals,
          FUN = "/"
        ) * 100
      )
    }
    
    stop(
      "Unknown abundance scale: ",
      scale,
      call. = FALSE
    )
  }
  
  if (normalise) {
    
    taxa <- normalize_to_percentage(
      taxa,
      taxa_scale
    )
    
    pathways <- normalize_to_percentage(
      pathways,
      pathways_scale
    )
  }
  
  ## ================================================================
  ## 9. QC report
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
  ## 10. Return MiCARA input object
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
 