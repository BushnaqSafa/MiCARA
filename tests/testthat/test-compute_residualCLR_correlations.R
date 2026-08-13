# 1. Helper to create lightweight, valid mock objects for tests
make_mock_cor_inputs <- function(n_samples = 20, n_taxa = 5, n_pathways = 5) {
  sample_ids <- paste0("Sample_", seq_len(n_samples))
  disease_labels <- rep(c("Control", "DiseaseA"), length.out = n_samples)

  # Metadata
  metadata <- data.frame(
    sample_id = sample_ids,
    disease = disease_labels,
    age = sample(20:60, n_samples, replace = TRUE),
    BMI = runif(n_samples, 18, 30),
    stringsAsFactors = FALSE
  )

  # Count matrices
  taxa_mat <- matrix(
    rpois(n_taxa * n_samples, lambda = 50),
    nrow = n_taxa, ncol = n_samples,
    dimnames = list(paste0("Taxon_", seq_len(n_taxa)), sample_ids)
  )
  pathway_mat <- matrix(
    rpois(n_pathways * n_samples, lambda = 100),
    nrow = n_pathways, ncol = n_samples,
    dimnames = list(paste0("Pathway_", seq_len(n_pathways)), sample_ids)
  )

  micara_obj <- structure(
    list(
      taxa = taxa_mat,
      pathways = pathway_mat,
      metadata = metadata,
      confounders = c("age", "BMI")
    ),
    class = "micara_input"
  )

  # Mock differential abundance objects
  taxa_diffab <- structure(
    list(
      feature_type = "taxa",
      significant = data.frame(
        feature = paste0("Taxon_", 1:3),
        disease = "DiseaseA",
        lfc = c(1.5, -0.8, 2.1),
        stringsAsFactors = FALSE
      )
    ),
    class = "micara_diffab"
  )

  pathway_diffab <- structure(
    list(
      feature_type = "pathways",
      significant = data.frame(
        feature = paste0("Pathway_", 1:3),
        disease = "DiseaseA",
        lfc = c(-1.2, 0.9, 1.1),
        stringsAsFactors = FALSE
      )
    ),
    class = "micara_diffab"
  )

  diffab_list <- structure(
    list(taxa = taxa_diffab, pathways = pathway_diffab),
    class = "micara_diffab_list"
  )

  list(
    micara_obj = micara_obj,
    taxa_diffab = taxa_diffab,
    pathway_diffab = pathway_diffab,
    diffab_list = diffab_list
  )
}

# 2. Input and class validation
test_that("compute_residualCLR_correlations validates input object classes", {
  inputs <- make_mock_cor_inputs()

  expect_error(
    compute_residualCLR_correlations(
      micara_obj = list(),
      diffab_obj = inputs$diffab_list,
      verbose = FALSE
    ),
    "'micara_obj' must inherit from 'micara_input'"
  )
})

# 3. Parameter range checks
test_that("compute_residualCLR_correlations checks numeric parameters and cutoffs", {
  inputs <- make_mock_cor_inputs()

  # q_cutoff bounds
  expect_error(
    compute_residualCLR_correlations(inputs$micara_obj, inputs$diffab_list, q_cutoff = 0, verbose = FALSE),
    "'q_cutoff' must be a single number in \\(0, 1\\]"
  )
  expect_error(
    compute_residualCLR_correlations(inputs$micara_obj, inputs$diffab_list, q_cutoff = 1.5, verbose = FALSE),
    "'q_cutoff' must be a single number in \\(0, 1\\]"
  )

  # r_cutoff bounds (allows 0)
  expect_error(
    compute_residualCLR_correlations(inputs$micara_obj, inputs$diffab_list, r_cutoff = -0.1, verbose = FALSE),
    "'r_cutoff' must be a single number in \\[0, 1\\]"
  )

  # Integer parameters
  expect_error(
    compute_residualCLR_correlations(inputs$micara_obj, inputs$diffab_list, min_samples_per_disease = 2, verbose = FALSE),
    "'min_samples_per_disease' must be one number >= 3"
  )
  expect_error(
    compute_residualCLR_correlations(inputs$micara_obj, inputs$diffab_list, min_residual_df = 0, verbose = FALSE),
    "'min_residual_df' must be one positive number"
  )
  expect_error(
    compute_residualCLR_correlations(inputs$micara_obj, inputs$diffab_list, retain_all_pairs = "invalid", verbose = FALSE),
    "'retain_all_pairs' must be TRUE or FALSE"
  )
})

# 4. Differential Abundance Object Resolution
test_that("compute_residualCLR_correlations resolves diffab input combinations correctly", {
  inputs <- make_mock_cor_inputs()

  # Cannot supply both diffab_obj and separate taxa/pathway objects
  expect_error(
    compute_residualCLR_correlations(
      inputs$micara_obj,
      diffab_obj = inputs$diffab_list,
      taxa_diffab = inputs$taxa_diffab,
      verbose = FALSE
    ),
    "Supply either a combined 'diffab_obj' or separate taxa/pathway objects, not both"
  )

  # Cannot supply single diffab object twice
  expect_error(
    compute_residualCLR_correlations(
      inputs$micara_obj,
      diffab_obj = inputs$taxa_diffab,
      taxa_diffab = inputs$taxa_diffab,
      verbose = FALSE
    ),
    "Taxa differential-abundance results were supplied twice"
  )

  # Invalid diffab class
  expect_error(
    compute_residualCLR_correlations(
      inputs$micara_obj,
      diffab_obj = list(),
      verbose = FALSE
    ),
    "'diffab_obj' must be a 'micara_diffab' or 'micara_diffab_list' object"
  )
})


# 5. Metadata & Covariate Integrity
test_that("compute_residualCLR_correlations validates metadata integrity", {
  inputs <- make_mock_cor_inputs()

  # Missing column
  bad_obj <- inputs$micara_obj
  bad_obj$metadata$disease <- NULL
  expect_error(
    compute_residualCLR_correlations(bad_obj, inputs$diffab_list, verbose = FALSE),
    "Missing required metadata column\\(s\\): disease"
  )

  # Missing or empty sample IDs
  bad_obj2 <- inputs$micara_obj
  bad_obj2$metadata$sample_id[1] <- NA
  expect_error(
    compute_residualCLR_correlations(bad_obj2, inputs$diffab_list, verbose = FALSE),
    "'metadata\\$sample_id' contains missing or empty values"
  )

  # Duplicate sample IDs
  bad_obj3 <- inputs$micara_obj
  bad_obj3$metadata$sample_id[2] <- bad_obj3$metadata$sample_id[1]
  expect_error(
    compute_residualCLR_correlations(bad_obj3, inputs$diffab_list, verbose = FALSE),
    "'metadata\\$sample_id' must contain unique sample IDs"
  )
})

# 6. CLR Transformation & Pseudocounts
test_that("compute_residualCLR_correlations checks covariate availability", {
  inputs <- make_mock_cor_inputs()

  expect_error(
    compute_residualCLR_correlations(
      inputs$micara_obj,
      inputs$diffab_list,
      covariates = c("non_existent_covariate"),
      verbose = FALSE
    ),
    "Covariate\\(s\\) not found in metadata: non_existent_covariate"
  )
})


test_that("compute_residualCLR_correlations parses pseudocount inputs accurately", {
  inputs <- make_mock_cor_inputs()

  # Single numeric pseudocount
  res1 <- compute_residualCLR_correlations(
    inputs$micara_obj,
    inputs$diffab_list,
    pseudocount = 0.5,
    verbose = FALSE
  )
  expect_equal(attr(res1, "settings")$taxa_pseudocount, 0.5)
  expect_equal(attr(res1, "settings")$pathway_pseudocount, 0.5)

  # Named list pseudocount
  res2 <- compute_residualCLR_correlations(
    inputs$micara_obj,
    inputs$diffab_list,
    pseudocount = list(taxa = 0.1, pathways = 0.2),
    verbose = FALSE
  )
  expect_equal(attr(res2, "settings")$taxa_pseudocount, 0.1)
  expect_equal(attr(res2, "settings")$pathway_pseudocount, 0.2)

  # Invalid pseudocount structure
  expect_error(
    compute_residualCLR_correlations(
      inputs$micara_obj,
      inputs$diffab_list,
      pseudocount = list(taxa = -1, pathways = 0.2),
      verbose = FALSE
    ),
    "Named pseudocount values must be NULL or single positive numbers"
  )
})

# 7. Disease alignment scoping
test_that("compute_residualCLR_correlations evaluates disease scope and warnings", {
  inputs <- make_mock_cor_inputs()

  # Warns when user requests a disease lacking matching features
  expect_warning(
    compute_residualCLR_correlations(
      inputs$micara_obj,
      inputs$diffab_list,
      diseases = c("DiseaseA", "NonExistentDisease"),
      verbose = FALSE
    ),
    "The following requested disease\\(s\\) do not have both significant taxa and significant pathways"
  )

  # Throws error if no diseases have both taxa and pathways
  empty_taxa_diffab <- inputs$taxa_diffab
  empty_taxa_diffab$significant <- empty_taxa_diffab$significant[0, ]

  expect_error(
    compute_residualCLR_correlations(
      inputs$micara_obj,
      taxa_diffab = empty_taxa_diffab,
      pathway_diffab = inputs$pathway_diffab,
      verbose = FALSE
    ),
    "No diseases have both significant taxa and significant pathways"
  )
})

# 8. End-to-End Interaction Computation
test_that("compute_residualCLR_correlations computes interactions successfully", {
  skip_if_not_installed("Hmisc")
  inputs <- make_mock_cor_inputs(n_samples = 30)

  res <- compute_residualCLR_correlations(
    inputs$micara_obj,
    diffab_obj = inputs$diffab_list,
    covariates = c("age", "BMI"),
    q_cutoff = 1.0, # Accept all correlations to test link generation
    r_cutoff = 0.0,
    retain_all_pairs = TRUE,
    verbose = FALSE
  )

  expect_s3_class(res, "micara_interactions")
  expect_type(res, "list")
  expect_named(res, "DiseaseA")

  disease_res <- res[["DiseaseA"]]
  expect_equal(disease_res$status, "analysed")
  expect_true(is.data.frame(disease_res$links))
  expect_true(is.data.frame(disease_res$all_links))
  expect_true(is.data.frame(disease_res$node_directions))

  # Check summary and settings attributes
  summary_tbl <- attr(res, "summary")
  settings_tbl <- attr(res, "settings")

  expect_true(is.data.frame(summary_tbl))
  expect_equal(summary_tbl$disease, "DiseaseA")
  expect_equal(settings_tbl$cor_method, "spearman")
})

# 9. Covariate Residualization & Degree-of-Freedom Guardrails
test_that("compute_residualCLR_correlations handles zero-variance covariates and degrees of freedom", {
  skip_if_not_installed("Hmisc")
  inputs <- make_mock_cor_inputs()

  # Add zero variance covariate
  inputs$micara_obj$metadata$constant_cov <- 10

  res <- compute_residualCLR_correlations(
    inputs$micara_obj,
    inputs$diffab_list,
    covariates = c("age", "constant_cov"),
    verbose = FALSE
  )

  disease_res <- res[["DiseaseA"]]
  expect_true(any(grepl("constant_cov", disease_res$covariates_dropped)))
  expect_true("age" %in% disease_res$covariates_used)

  # Force residual DF failure
  expect_message(
    res_df <- compute_residualCLR_correlations(
      inputs$micara_obj,
      inputs$diffab_list,
      min_residual_df = 100, # Impossibly high DF requirement
      verbose = TRUE
    ),
    "skipped"
  )
  expect_equal(res_df[["DiseaseA"]]$status, "skipped")
})

# 10. Optional pair retension
test_that("compute_residualCLR_correlations respects retain_all_pairs = FALSE", {
  skip_if_not_installed("Hmisc")
  inputs <- make_mock_cor_inputs()

  res <- compute_residualCLR_correlations(
    inputs$micara_obj,
    inputs$diffab_list,
    retain_all_pairs = FALSE,
    verbose = FALSE
  )

  expect_null(res[["DiseaseA"]]$all_links)
})
