#' #' Validate and Preprocess MiCARA Input Tables
#'
#' Performs quality control on taxonomic and functional pathway abundance
#' tables and accompanying sample metadata. Checks sample ID alignment,
#' detects abundance scales (counts vs. relative abundances), back-calculates 
#' estimated raw integer counts when necessary using sequencing depth, 
#' filters rare disease groups, and evaluates metadata covariate missingness.
#' Automatically detects and drops leading serial or index columns 
#' (e.g., \code{1:N}, \code{X}, \code{Unnamed: 0}) if present.
#'
#' @param taxa data.frame of taxon abundances with samples as columns and 
#'   taxa as rows. Rownames can be set directly or provided in a leading 
#'   non-numeric column. Leading serial or index columns (e.g., \code{1:N}, \code{X}) 
#'   are automatically detected and dropped.
#' @param pathways data.frame of pathway abundances with samples as columns 
#'   and pathways as rows. Follows the same structure as \code{taxa}.
#' @param metadata data.frame with one row per sample. Must contain 
#'   \code{sample_id} and \code{disease}. Optional recognized columns include 
#'   sequencing depth (\code{number_reads}, \code{sequencing_depth}, 
#'   \code{total_reads}, \code{n_reads}, or \code{depth}—used to reconstruct 
#'   raw counts if relative abundance is provided), \code{study_name}, 
#'   \code{age}, \code{age_category}, \code{gender}, and \code{BMI}.
#' @param disease_min_samples Minimum number of samples required to retain 
#'   a disease group (default 10).
#' @param metadata_missing_cutoff Maximum allowed fraction of missing 
#'   values for an optional covariate to be retained (default 0.30).
#' @param verbose Logical; if \code{TRUE} (default), prints a summary report 
#'   to the console.
#' @param disease_col Character string specifying the metadata column containing disease cohort labels.
#' @param normalise Logical indicating whether feature normalization should be performed. Default is TRUE.
#' 
#' 
#' 
#' 
#' @return An object of class \code{micara_input}: a list containing:
#'   \item{taxa}{Processed data frame of taxon abundances (reconstructed to 
#'     estimated raw counts if relative abundances were supplied alongside 
#'     a sequencing depth column in metadata).}
#'   \item{pathways}{Processed data frame of pathway abundances.}
#'   \item{metadata}{Filtered and aligned metadata data frame.}
#'   \item{abundance_scale}{A list containing \code{detected} (initial scale) 
#'     and \code{final} (resulting scale after processing).}
#'   \item{confounders}{Character vector of retained covariate names.}
#'   \item{removed_covariates}{Character vector of covariates excluded due to missingness.}
#'   \item{removed_diseases}{Character vector of disease groups excluded due to low sample size.}
#'
#' @export  


  validate_inputs <- function(
    taxa,
    pathways,
    metadata,
    disease_col = "disease",
    disease_min_samples     = 10,
    metadata_missing_cutoff = 0.30,
    normalise               = TRUE,
    verbose                 = TRUE
  ) {

    ## ------------------------------------------------------------------
    ## 1a. Object types & Auto-Coercion
    ## ------------------------------------------------------------------
    
    coerce_to_df <- function(x, name) {
      if (is.matrix(x) || is.data.frame(x) || inherits(x, c("DFrame", "DataFrame", "tbl_df", "tbl", "data.table"))) {
        return(as.data.frame(x))
      }
      stop(sprintf("'%s' must be a matrix, data.frame, tibble, or Bioconductor DFrame (received class: %s)", name, class(x)[1]))
    }
    
    taxa     <- coerce_to_df(taxa, "taxa")
    pathways <- coerce_to_df(pathways, "pathways")
    metadata <- coerce_to_df(metadata, "metadata")
    
    ## ------------------------------------------------------------------
    ## 1b. Helper: Drop serial numbers & handle feature row names
    ## ------------------------------------------------------------------
    move_first_col_to_rownames <- function(df, label) {
      df <- as.data.frame(df, stringsAsFactors = FALSE)
      if (ncol(df) == 0) return(df)
      
      # Step A: Detect and drop serial / index column from Column 1
      c1     <- df[[1]]
      cn     <- colnames(df)[1]
      is_seq <- is.numeric(c1) && (identical(as.integer(c1), seq_len(nrow(df))) || identical(as.integer(c1), 0L:(nrow(df) - 1L)))
      is_idx <- grepl("^(X|Unnamed|index|serial|row_?num)", cn, ignore.case = TRUE)
      
      if (is_seq || is_idx) {
        message(sprintf("[MiCARA] Automatically dropped serial/index column '%s' from %s.", cn, label))
        df <- df[, -1, drop = FALSE]
      }
      
      # Step B: Evaluate remaining columns
      is_numeric_col <- vapply(df, is.numeric, logical(1))
      
      # All numeric -> already clean matrix with rownames attached
      if (all(is_numeric_col)) {
        return(df)
      }
      
      # Non-numeric 1st column + numeric rest -> move feature names to rownames
      if (!is_numeric_col[1] && all(is_numeric_col[-1])) {
        feature_names <- as.character(df[[1]])
        
        if (any(duplicated(feature_names))) {
          stop(
            "'", label, "' first column contains duplicate names; ",
            "cannot use as unique row identifiers.",
            call. = FALSE
          )
        }
        
        df <- df[, -1, drop = FALSE]
        rownames(df) <- feature_names
        return(df)
      }
      
      stop(
        "'", label, "' contains non-numeric columns that could not be auto-resolved. ",
        "Provide either (a) rownames = feature names with all columns numeric, or ",
        "(b) a single leading non-numeric column of feature names followed by numeric sample columns.",
        call. = FALSE
      )
    }
    
    # Clean taxa and pathways
    taxa     <- move_first_col_to_rownames(taxa,     "taxa")
    pathways <- move_first_col_to_rownames(pathways, "pathways")
    
    ## ------------------------------------------------------------------
    ## 1c. Disease / Group column validation & standardization
    ## ------------------------------------------------------------------
    
    if (!disease_col %in% colnames(metadata)) {
      stop(sprintf(
        "Specified `disease_col` ('%s') was not found in metadata columns. Available columns: %s",
        disease_col,
        paste(sQuote(colnames(metadata)), collapse = ", ")
      ))
    }
    
    # Standardize internally so rest of function can always rely on metadata$disease
    metadata$disease <- metadata[[disease_col]]
    
    
    ## ------------------------------------------------------------------
    ## 2. Flexible feature-name handling + numeric content check
    ## cleans up imported spreadsheets/data frames so that features (like taxa names or pathway IDs) become row names, leaving only numeric sample measurements in the matrix.
    ##
    ## Accepts either:
    ##   (a) rownames already set to taxon/pathway names, all columns numeric
    ##   (b) first column holds taxon/pathway names as a data column, with
    ##       every other column numeric (typical of a TSV read via
    ##       read.delim without row.names = 1)
    ## Any other non-numeric column configuration is rejected.
    ## ------------------------------------------------------------------
    
    move_first_col_to_rownames <- function(df, label) {
      
      is_numeric_col <- vapply(df, is.numeric, logical(1))
      
      if (all(is_numeric_col)) {
        return(df)
      }
      
      if (!is_numeric_col[1] && all(is_numeric_col[-1])) {
        
        feature_names <- as.character(df[[1]])
        
        if (any(duplicated(feature_names))) {
          stop(
            "'", label, "' first column contains duplicate names; ",
            "cannot use as unique row identifiers."
          )
        }
        
        df <- df[, -1, drop = FALSE]
        rownames(df) <- feature_names
        
        return(df)
      }
      
      stop(
        "'", label, "' contains non-numeric columns that could not be ",
        "auto-resolved. Provide either (a) rownames = feature names with ",
        "all columns numeric, or (b) a single leading non-numeric column ",
        "of feature names followed by numeric sample columns."
      )
    }
    
    taxa     <- move_first_col_to_rownames(taxa,     "taxa")
    pathways <- move_first_col_to_rownames(pathways, "pathways")
    
    ## ------------------------------------------------------------------
    ## 2b. Auto-cleaning MetaPhlAn taxa & HUMAnN pathways
    ## ------------------------------------------------------------------
    
    # Auto-clean MetaPhlAn taxonomic names if detected
    if (any(grepl("s__", rownames(taxa)))) {
      if (verbose) {
        message("[Info] MetaPhlAn taxonomic prefixes detected. Extracting species and cleaning names...")
      }
      keep_taxa <- grepl("s__", rownames(taxa)) & !grepl("t__", rownames(taxa))
      taxa <- taxa[keep_taxa, , drop = FALSE]
      rownames(taxa) <- chartr("_", " ", sub(".*s__", "", rownames(taxa)))
    }
    
    # Auto-clean HUMAnN pathway rows if detected
    if (any(grepl("\\|", rownames(pathways))) || any(grepl("UNINTEGRATED|UNMAPPED", rownames(pathways)))) {
      if (verbose) {
        message("[Info] HUMAnN stratified/unmapped pathways detected. Extracting community-level pathways...")
      }
      keep_pwy <- !grepl("\\|", rownames(pathways)) & !grepl("UNINTEGRATED|UNMAPPED", rownames(pathways))
      pathways <- pathways[keep_pwy, , drop = FALSE]
    }
    
    ## ------------------------------------------------------------------
    ## 3. Sample ID alignment (taxa vs pathways)
    ## ------------------------------------------------------------------
    
    if (!identical(colnames(taxa), colnames(pathways))) {
      stop("Sample IDs in 'taxa' and 'pathways' are not identical and in the same order.")
    }
    
    if (!all(colnames(taxa) %in% metadata$sample_id)) {
      missing_ids <- setdiff(colnames(taxa), metadata$sample_id)
      stop(
        length(missing_ids), " sample ID(s) present in taxa/pathways ",
        "but absent from metadata (e.g. ", paste(utils::head(missing_ids, 3), collapse = ", "), ")"
      )
    }
    
    ## reorder metadata to match taxa/pathways column order
    metadata <- metadata[match(colnames(taxa), metadata$sample_id), ]
    
    ## defensive re-check after reordering
    if (!identical(colnames(taxa), metadata$sample_id)) {
      stop("Internal alignment error: metadata$sample_id does not match taxa columns after reordering.")
    }
    
    ## ------------------------------------------------------------------
    ## 4. Mandatory metadata columns
    ## ------------------------------------------------------------------
    
    mandatory <- c("sample_id", "disease")
    missing_mandatory <- setdiff(mandatory, names(metadata))
    
    if (length(missing_mandatory) > 0) {
      stop("Missing mandatory metadata column(s): ", paste(missing_mandatory, collapse = ", "))
    }
    
    # ------------------------------------------------------------------
    # 4b. Flexible detection for total read counts / depth column
    # ------------------------------------------------------------------
    read_count_aliases <- c(
      "number_reads", "sequencing_depth", "total_reads", 
      "read_count", "read_counts", "n_reads", "depth", "reads"
    )
    
    # Search for matching column names (case-insensitive)
    metadata_lower <- tolower(names(metadata))
    matched_alias_idx <- which(metadata_lower %in% read_count_aliases)
    
    if (length(matched_alias_idx) > 0) {
      found_col <- names(metadata)[matched_alias_idx[1]]
      
      # Rename to standardized "number_reads" internal key
      if (found_col != "number_reads") {
        names(metadata)[matched_alias_idx[1]] <- "number_reads"
        if (verbose) {
          message("[Info] Detected sequencing depth column '", found_col, "' and mapped it to 'number_reads'.")
        }
      }
    } else if (verbose) {
      message(
        "[Info] No total read count column detected in metadata. ",
        "If you have relative abundance tables and want ANCOM-BC2 raw count reconstruction, ",
        "include a column named 'number_reads' or 'sequencing_depth'."
      )
    }
    
    ## ------------------------------------------------------------------
    ## 5. Disease filtering 
    ## ------------------------------------------------------------------
    
    disease_counts  <- table(metadata$disease)
    keep_disease    <- names(disease_counts[disease_counts >= disease_min_samples])
    removed_disease <- setdiff(names(disease_counts), keep_disease)
    
    keep_samples <- metadata$disease %in% keep_disease
    
    n_removed_samples <- sum(!keep_samples)
    
    metadata <- metadata[keep_samples, ]
    taxa     <- taxa[,     keep_samples, drop = FALSE]
    pathways <- pathways[, keep_samples, drop = FALSE]
    
    ## ------------------------------------------------------------------
    ## 6. Optional covariate checking and missing-data filtering
    ## ------------------------------------------------------------------
    
    expected_optional <- c("study_name", "age", "age_category", "gender", "BMI")
    # 1. Map lower-case expected names to their standard casing
    expected_map <- stats::setNames(expected_optional, tolower(expected_optional))
    # 2. Find matching columns in metadata (case-insensitive)
    matched_idx  <- match(tolower(names(metadata)), names(expected_map))
    # 3. Rename any matched columns in metadata to standard casing
    names(metadata)[!is.na(matched_idx)] <- expected_map[matched_idx[!is.na(matched_idx)]]
    # 4. Filter for available expected covariates (now matching exact casing)
    available_optional <- intersect(expected_optional, names(metadata))
    # 5. Check missingness on available covariates
    missing_fraction <- vapply(
      available_optional,
      function(x) mean(is.na(metadata[[x]])),
      numeric(1)
    )
    
    keep_covariates   <- names(missing_fraction[missing_fraction <= metadata_missing_cutoff])
    remove_covariates <- names(missing_fraction[missing_fraction >  metadata_missing_cutoff])
    
  
    ## ------------------------------------------------------------------
    ## 7. HELPER: Detect abundance scale
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
    ## 8. HELPER: Back-calculate raw counts from relative abundance 
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
    ## 9. VERIFICATION & RECONSTRUCTION PIPELINE
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
    
    ## ------------------------------------------------------------------
    ## 10. Report
    ## ------------------------------------------------------------------
    
    if (verbose) {
      # Helper to construct clear scale reporting status
      get_scale_status <- function(initial_scale) {
        if (initial_scale == "counts") {
          return("Raw counts (Unmodified)")
        } else if (initial_scale %in% c("percentage", "proportion")) {
          if (has_reads) {
            return(paste0(initial_scale, " -> Reconstructed to estimated raw counts (via number_reads)"))
          } else {
            return(paste0(initial_scale, " -> Left unmodified (number_reads unavailable)"))
          }
        } else {
          return(paste0(initial_scale, " (Unmodified)"))
        }
      }
      
      cat("\n========== MiCARA Input Report ==========\n")
      cat("Samples retained:      ", ncol(taxa), "\n")
      cat("Samples removed:       ", n_removed_samples, "\n")
      cat("Diseases retained:     ", length(keep_disease), "\n")
      if (length(removed_disease) > 0) {
        cat("Diseases removed:      ", paste(removed_disease, collapse = ", "), 
            paste0("(n < ", disease_min_samples, ")\n"))
      } else {
        cat("Diseases removed:      none\n")
      }
      
      cat("Taxa abundance scale:  ", get_scale_status(taxa_scale), "\n")
      cat("Pathway scale:         ", get_scale_status(pathways_scale), "\n")
      
      cat("Covariates retained:   ", 
          if (length(keep_covariates) > 0) paste(keep_covariates, collapse = ", ") else "none", "\n")
      cat("Covariates removed:    ", 
          if (length(remove_covariates) > 0) {
            paste0(paste(remove_covariates, collapse = ", "), " (>", metadata_missing_cutoff * 100, "% missing)")
          } else {
            "none"
          }, "\n")
      cat("------------------------------------------\n")
    }
    
    ## ------------------------------------------------------------------
    ## 11. Return
    ## ------------------------------------------------------------------
    # Determine final output scale based on whether reconstruction occurred
    taxa_final_scale     <- if (taxa_scale %in% c("percentage", "proportion") && has_reads) "counts" else taxa_scale
    pathways_final_scale <- if (pathways_scale %in% c("percentage", "proportion") && has_reads) "counts" else pathways_scale
    
    # Final safety check: ensure metadata rows match taxa column order exactly
    metadata <- metadata[match(colnames(taxa), metadata$sample_id), ]
    
    structure(
      list(
        taxa               = taxa,
        pathways           = pathways,
        metadata           = metadata,
        abundance_type    = list(
          detected = list(taxa = taxa_scale, pathways = pathways_scale),
          final    = list(taxa = taxa_final_scale, pathways = pathways_final_scale)
        ),
        confounders        = if (is.null(keep_covariates)) character(0) else keep_covariates,
        removed_covariates = if (is.null(remove_covariates)) character(0) else remove_covariates,
        removed_diseases   = if (is.null(removed_disease)) character(0) else removed_disease
      ),
      class = "micara_input"
    )
  }
  

 