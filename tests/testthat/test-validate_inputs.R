## tests/testthat/test-validate_inputs.R
##
## Synthetic datasets are constructed here rather than loaded from disk so
## that tests are self-contained, fast, and do not depend on external files.



## ------------------------------------------------------------------------
## 1. Happy-path test (retain percentage input without read counts)
## ------------------------------------------------------------------------

test_that("valid percentage-scale input is retained when read counts are unavailable", {
  
  d <- make_test_data(abundance_scale = "percentage")
  
  result <- suppressWarnings(
    validate_inputs(
      d$taxa, d$pathways, d$metadata,
      verbose = FALSE
    )
  )
  
  expect_s3_class(result, "micara_input")
  
  expect_equal(result$abundance_type$detected$taxa, "percentage")
  expect_equal(result$abundance_type$detected$pathways, "percentage")
  
  expect_equal(result$abundance_type$final$taxa, "percentage")
  expect_equal(result$abundance_type$final$pathways, "percentage")
  
  expect_equal(result$taxa, d$taxa)
  expect_equal(result$pathways, d$pathways)
  expect_equal(result$metadata, d$metadata)
  
  expect_true(all(abs(colSums(result$taxa) - 100) < 1e-6))
  expect_true(all(abs(colSums(result$pathways) - 100) < 1e-6))
})



test_that("relative abundance without number_reads issues a warning", {
  
  d <- make_test_data(abundance_scale = "percentage")
  
  suppressWarnings(expect_warning(
    validate_inputs(
      d$taxa,
      d$pathways,
      d$metadata,
      verbose = FALSE),
    "Input is in 'percentage' scale")
  )
})

## ------------------------------------------------------------------------
## 1. Test counts reconstruction when number_reads is available
## ------------------------------------------------------------------------

#Percentage reconstruction test
test_that("percentage-scale input is reconstructed to estimated counts when read counts are available", {
  
  d <- make_test_data(
    abundance_scale = "percentage",
    add_number_reads = TRUE
  )
  
  expect_message(
    suppressWarnings(result <- validate_inputs(
      d$taxa,
      d$pathways,
      d$metadata,
      verbose = FALSE
    )),
    "Back-calculated estimated counts"
  )
  
  expect_equal(result$abundance_type$detected$taxa, "percentage")
  expect_equal(result$abundance_type$detected$pathways, "percentage")
  
  expect_equal(result$abundance_type$final$taxa, "counts")
  expect_equal(result$abundance_type$final$pathways, "counts")
  
  expect_equal(
    result$taxa,
    expected_reconstructed_counts(
      d$taxa,
      "percentage",
      d$metadata$number_reads
    )
  )
  
  expect_equal(
    result$pathways,
    expected_reconstructed_counts(
      d$pathways,
      "percentage",
      d$metadata$number_reads
    )
  )
})

#Proportion reconstruction test
test_that("proportion-scale input is reconstructed to estimated counts when read counts are available", {
  
  d <- make_test_data(
    abundance_scale = "proportion",
    add_number_reads = TRUE
  )
  
  expect_message(suppressWarnings(result <- validate_inputs(
    d$taxa,
    d$pathways,
    d$metadata,
    verbose = FALSE))
    ,
    "Back-calculated estimated counts"
  )
  
  expect_equal(result$abundance_type$detected$taxa, "proportion")
  expect_equal(result$abundance_type$detected$pathways, "proportion")
  
  expect_equal(result$abundance_type$final$taxa, "counts")
  expect_equal(result$abundance_type$final$pathways, "counts")
  
  expect_equal(
    result$taxa,
    expected_reconstructed_counts(
      d$taxa,
      "proportion",
      d$metadata$number_reads
    )
  )
  
  expect_equal(
    result$pathways,
    expected_reconstructed_counts(
      d$pathways,
      "proportion",
      d$metadata$number_reads
    )
  )
})

## ------------------------------------------------------------------------
## 2. Object type check
## ------------------------------------------------------------------------

test_that("non-data.frame inputs are rejected", {
  
  d <- make_test_data()
  
  expect_error(
    validate_inputs(as.matrix(d$taxa), d$pathways, d$metadata, verbose = FALSE),
    "must be a data.frame"
  )
  
  expect_error(
    validate_inputs(d$taxa, as.matrix(d$pathways), d$metadata, verbose = FALSE),
    "must be a data.frame"
  )
  
  expect_error(
    validate_inputs(d$taxa, d$pathways, as.matrix(d$metadata), verbose = FALSE),
    "must be a data.frame"
  )
})


## ------------------------------------------------------------------------
## 3. Non-numeric content
## ------------------------------------------------------------------------

test_that("non-numeric columns in taxa/pathways are caught", {
  
  d <- make_test_data()
  d$taxa$extra_text_col <- "not_a_number"
  
  expect_error(
    validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE),
    "contains non-numeric columns"
  )
})


## ------------------------------------------------------------------------
## 4. Sample ID misalignment between taxa and pathways
## ------------------------------------------------------------------------

test_that("misaligned taxa/pathways columns are detected", {
  
  d <- make_test_data(misalign_pathways = TRUE)
  
  expect_error(
    validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE),
    "not identical and in the same order"
  )
})


## ------------------------------------------------------------------------
## 5. Sample present in taxa but missing from metadata
## ------------------------------------------------------------------------

test_that("sample missing from metadata triggers informative error", {
  
  d <- make_test_data(drop_metadata_id = "healthy_1")
  
  expect_error(
    validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE),
    "sample ID\\(s\\) present in taxa/pathways"
  )
})


## ------------------------------------------------------------------------
## 6. Mandatory metadata columns missing
## ------------------------------------------------------------------------

test_that("missing mandatory metadata column ('disease') is caught", {
  
  d <- make_test_data()
  d$metadata$disease <- NULL
  
  expect_error(
    validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE),
    "not found in metadata columns"
  )
})


## ------------------------------------------------------------------------
## 7. Disease filtering (n < disease_min_samples)
## ------------------------------------------------------------------------

test_that("diseases below the sample size threshold are removed", {
  d <- make_test_data(n_samples_per_disease = c(healthy = 15, T2D = 15, rare_disease = 3))
  
  suppressWarnings({
    result <- validate_inputs(d$taxa, d$pathways, d$metadata, disease_min_samples = 10, verbose = FALSE)
  })
  
  expect_false("rare_disease" %in% result$metadata$disease)
  expect_equal(result$removed_diseases, "rare_disease")
  expect_equal(ncol(result$taxa), 30)
})

## ------------------------------------------------------------------------
## 8. Covariate missingness cutoff
## ------------------------------------------------------------------------

test_that("covariates exceeding the missingness cutoff are dropped", {
  d <- make_test_data(covariate_missing_frac = 0.5)
  suppressWarnings({
    result <- validate_inputs(d$taxa, d$pathways, d$metadata, metadata_missing_cutoff = 0.30, verbose = FALSE)
  })
  expect_false("BMI" %in% result$confounders)
  expect_true("BMI" %in% result$removed_covariates)
})

test_that("covariates below the missingness cutoff are retained", {
  d <- make_test_data(covariate_missing_frac = 0.10)
  suppressWarnings({
    result <- validate_inputs(d$taxa, d$pathways, d$metadata, metadata_missing_cutoff = 0.30, verbose = FALSE)
  })
  expect_true("BMI" %in% result$confounders)
})

## ------------------------------------------------------------------------
## 9. Covariate near-miss naming (case mismatch)
## ------------------------------------------------------------------------

test_that("alternative cased covariate name in metadata is correctly detected", {
  d <- make_test_data(rename_bmi = TRUE)
  suppressWarnings({
    result <- validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
  })
  expect_true("bmi" %in% result$confounders || "BMI" %in% result$confounders)
})

## ------------------------------------------------------------------------
## 10. Abundance-scale detection and processing
## ------------------------------------------------------------------------

# Proportions test
test_that("proportion-scale input is correctly detected and left unchanged", {
  
  d <- make_test_data(abundance_scale = "proportion")
  suppressWarnings({
    result <- validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
  })
  
  expect_equal(result$abundance_type$detected$taxa, "proportion")
  expect_equal(result$abundance_type$detected$pathways, "proportion")
  
  ## No number_reads column is present, so reconstruction does not occur.
  expect_equal(result$abundance_type$final$taxa, "proportion")
  expect_equal(result$abundance_type$final$pathways, "proportion")
  
  expect_equal(result$taxa, d$taxa)
  expect_equal(result$pathways, d$pathways)
  
  expect_true(all(abs(colSums(result$taxa) - 1) < 1e-6))
  expect_true(all(abs(colSums(result$pathways) - 1) < 1e-6))
})



#Counts test
test_that("count-scale input is correctly detected and left unchanged", {
  
  d <- make_test_data(abundance_scale = "counts")
  
  result <- validate_inputs(
    d$taxa,
    d$pathways,
    d$metadata,
    verbose = FALSE
  )
  
  expect_equal(result$abundance_type$detected$taxa, "counts")
  expect_equal(result$abundance_type$detected$pathways, "counts")
  
  expect_equal(result$abundance_type$final$taxa, "counts")
  expect_equal(result$abundance_type$final$pathways, "counts")
  
  expect_equal(result$taxa, d$taxa)
  expect_equal(result$pathways, d$pathways)
})

#unrecognisable abundance scale
test_that("unrecognisable abundance scale raises an informative error", {
  
  d <- make_test_data(abundance_scale = "proportion")
  
  ## Produce column totals of 0.5:
  ## not percentage (~100), proportion (~1), or counts (>1).
  d$taxa <- d$taxa * 0.5
  
  expect_error(
    suppressWarnings(
      validate_inputs(
        d$taxa,
        d$pathways,
        d$metadata,
        verbose = FALSE
      )
    ),
    "Unable to determine abundance scale"
  )
})


## ------------------------------------------------------------------------
## 11. Return object structure
## ------------------------------------------------------------------------

test_that("returned object has the expected structure and class", {
  
  d <- make_test_data()
  suppressWarnings({result <- validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
  })
  
  expect_s3_class(result, "micara_input")
  expect_named(
    result, c("taxa", "pathways", "metadata", "abundance_type", "confounders", "removed_covariates", "removed_diseases")
  )
  expect_identical(colnames(result$taxa), colnames(result$pathways))
  expect_identical(colnames(result$taxa), result$metadata$sample_id)
})

## ------------------------------------------------------------------------
## 12. Edge case: Custom disease column standardization
## ------------------------------------------------------------------------

test_that("custom disease column name is accepted and standardized to disease", {
  d <- make_test_data()
  d$metadata$Group <- d$metadata$disease
  d$metadata$disease <- NULL
  
  result <- suppressWarnings(
    validate_inputs(
      d$taxa, d$pathways, d$metadata,
      disease_col = "Group",
      verbose = FALSE
    )
  )
  
  expect_true("disease" %in% colnames(result$metadata))
  expect_equal(result$metadata$disease, d$metadata$Group)
})


## ------------------------------------------------------------------------
## 13. First-column feature names handling
## ------------------------------------------------------------------------

test_that("feature names in the first column are converted to rownames", {
  d <- make_test_data()
  taxa_df <- cbind(feature_id = rownames(d$taxa), d$taxa)
  rownames(taxa_df) <- NULL
  
  result <- suppressWarnings(
    validate_inputs(taxa_df, d$pathways, d$metadata, verbose = FALSE)
  )
  
  expect_equal(rownames(result$taxa), paste0("taxon_", 1:10))
  expect_true(is.numeric(as.matrix(result$taxa)))
})

test_that("duplicate feature names in the first column raise an error", {
  d <- make_test_data()
  taxa_df <- cbind(feature_id = rownames(d$taxa), d$taxa)
  taxa_df$feature_id[2] <- taxa_df$feature_id[1]
  
  expect_error(
    validate_inputs(taxa_df, d$pathways, d$metadata, verbose = FALSE),
    "duplicate|unique"
  )
})


## ------------------------------------------------------------------------
## 14. Sequencing depth alias auto-detection
## ------------------------------------------------------------------------

test_that("sequencing depth column aliases are mapped to number_reads", {
  d <- make_test_data()
  d$metadata$read_count <- sample(10000:50000, nrow(d$metadata), replace = TRUE)
  
  result <- suppressWarnings(
    validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
  )
  
  expect_true("number_reads" %in% colnames(result$metadata))
  expect_equal(result$metadata$number_reads, d$metadata$read_count)
})









