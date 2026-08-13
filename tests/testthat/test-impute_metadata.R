test_that("impute_metadata validates S3 input class", {
  d <- make_test_data()

  # Must fail if passed a raw data frame or list rather than micara_input object
  expect_error(
    impute_metadata(d$metadata),
    "'micara_obj' must be the output of validate_inputs"
  )
})

test_that("impute_metadata exits early when no variables contain missing values", {
  d <- make_test_data()
  validated <- suppressWarnings(
    validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
  )
  if (is.null(validated$confounders)) validated$confounders <- c("age", "gender", "BMI")

  res <- impute_metadata(validated, verbose = FALSE)

  expect_s3_class(res, "micara_imputed")
  expect_false(res$imputation_performed)
  expect_null(res$mice_fit)
  expect_equal(res$metadata, validated$metadata)
})

test_that("impute_metadata blocks age_category as an imputation target", {
  d <- make_test_data()
  validated <- suppressWarnings(
    validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
  )

  expect_error(
    impute_metadata(validated, impute_vars = "age_category", verbose = FALSE),
    "'age_category' should not be used as an imputation target"
  )
})

test_that("impute_metadata throws an error if specified impute_vars do not exist", {
  d <- make_test_data()
  validated <- suppressWarnings(
    validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
  )

  expect_error(
    impute_metadata(validated, impute_vars = "non_existent_col", verbose = FALSE),
    "impute_vars not found in metadata"
  )
})

test_that("impute_metadata throws an error if disease contains missing values", {
  d <- make_test_data(covariate_missing_frac = 0.2)


  validated <- suppressWarnings(
    validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
  )
  if (is.null(validated$confounders)) validated$confounders <- c("age", "gender", "BMI")

  # Inject NA into disease after validation
  validated$metadata$disease[1] <- NA

  expect_error(
    impute_metadata(validated, impute_vars = "BMI", verbose = FALSE),
    "'disease' contains missing values"
  )
})

test_that("impute_metadata executes MICE successfully and fills missing values", {
  d <- make_test_data(covariate_missing_frac = 0.2)
  expect_true(any(is.na(d$metadata$BMI)))

  validated <- suppressWarnings(
    validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
  )
  if (is.null(validated$confounders)) validated$confounders <- c("age", "gender", "BMI")

  res <- impute_metadata(validated, impute_vars = "BMI", m = 2, verbose = FALSE)

  expect_s3_class(res, "micara_imputed")
  expect_true(res$imputation_performed)
  expect_false(any(is.na(res$metadata$BMI)))
  expect_s3_class(res$mice_fit, "mids")
})
