## ------------------------------------------------------------------------
## Helper to build a minimal valid dataset, with arguments to deliberately
## break specific aspects for individual tests.
## ------------------------------------------------------------------------

make_test_data <- function(
    n_samples_per_disease  = c(healthy = 15, T2D = 15),
    abundance_scale        = "percentage",
    add_number_reads       = FALSE,
    misalign_pathways      = FALSE,
    drop_metadata_id       = NULL,
    covariate_missing_frac = 0,
    rename_bmi             = FALSE
) {
  
  sample_ids <- unlist(lapply(names(n_samples_per_disease), function(d) {
    paste0(d, "_", seq_len(n_samples_per_disease[[d]]))
  }))
  
  disease_vec <- rep(names(n_samples_per_disease), times = n_samples_per_disease)
  
  n <- length(sample_ids)
  
  set.seed(42)
  
  ## --- taxa table: 10 taxa x n samples ---
  raw <- matrix(runif(10 * n, 0, 10), nrow = 10)
  
  if (abundance_scale == "percentage") {
    taxa <- sweep(raw, 2, colSums(raw), "/") * 100
  } else if (abundance_scale == "proportion") {
    taxa <- sweep(raw, 2, colSums(raw), "/")
  } else if (abundance_scale == "counts") {
    taxa <- round(raw * 10000)
  }
  
  taxa <- as.data.frame(taxa)
  colnames(taxa) <- sample_ids
  rownames(taxa) <- paste0("taxon_", seq_len(10))
  
  ## --- pathways table: 6 pathways x n samples (same scale as taxa) ---
  raw_p <- matrix(runif(6 * n, 0, 10), nrow = 6)
  
  if (abundance_scale == "percentage") {
    pathways <- sweep(raw_p, 2, colSums(raw_p), "/") * 100
  } else if (abundance_scale == "proportion") {
    pathways <- sweep(raw_p, 2, colSums(raw_p), "/")
  } else if (abundance_scale == "counts") {
    pathways <- round(raw_p * 10000)
  }
  
  pathways <- as.data.frame(pathways)
  colnames(pathways) <- sample_ids
  rownames(pathways) <- paste0("pathway_", seq_len(6))
  
  if (misalign_pathways) {
    ## shuffle pathway columns so order no longer matches taxa
    pathways <- pathways[, rev(seq_len(n))]
  }
  
  ## --- metadata ---
  metadata <- data.frame(
    sample_id    = sample_ids,
    disease      = disease_vec,
    study_name   = sample(c("studyA", "studyB"), n, replace = TRUE),
    age          = sample(20:80, n, replace = TRUE),
    age_category = sample(c("adult", "senior"), n, replace = TRUE),
    gender       = sample(c("male", "female"), n, replace = TRUE),
    BMI          = round(rnorm(n, 25, 4), 1),
    stringsAsFactors = FALSE
  )
  
  if (add_number_reads) {
    metadata$number_reads <- sample(
      c(5000, 10000, 15000, 20000),
      n,
      replace = TRUE
    )
  }
  
  if (covariate_missing_frac > 0) {
    n_na <- floor(covariate_missing_frac * n)
    metadata$BMI[sample(seq_len(n), n_na)] <- NA
  }
  
  if (rename_bmi) {
    names(metadata)[names(metadata) == "BMI"] <- "bmi"
  }
  
  if (!is.null(drop_metadata_id)) {
    metadata <- metadata[metadata$sample_id != drop_metadata_id, ]
  }
  
  list(taxa = taxa, pathways = pathways, metadata = metadata)
}


## Helper to compute the expected reconstructed counts.
expected_reconstructed_counts <- function(mat, scale, sample_reads) {
  
  if (scale == "counts") {
    return(mat)
  }
  
  if (scale == "percentage") {
    prop_mat <- mat / 100
  } else if (scale == "proportion") {
    prop_mat <- mat
  } else {
    stop(
      "scale must be one of 'counts', 'percentage', or 'proportion'.",
      call. = FALSE
    )
  }
  
  round(sweep(prop_mat, 2, sample_reads, FUN = "*"))
}