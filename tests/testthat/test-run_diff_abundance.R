skip_if_not_installed("microbiome")
skip_if_not_installed("ANCOMBC")


test_that("run_diff_abundance validates micara_obj class", {
    expect_error(
        run_diff_abundance(list(taxa = data.frame())),
        "'micara_obj' must be an object of class 'micara_input'"
    )
})

test_that("run_diff_abundance validates feature_type parameter", {
    d <- make_test_data(abundance_scale = "counts")
    validated <- suppressWarnings(
        validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
    )

    expect_error(
        run_diff_abundance(validated, feature_type = "invalid_feature", verbose = FALSE),
        "Invalid feature_type specified"
    )
})

test_that("run_diff_abundance throws error when reference_level is absent in disease", {
    d <- make_test_data(abundance_scale = "counts")
    validated <- suppressWarnings(
        validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
    )

    expect_error(
        run_diff_abundance(validated, feature_type = "taxa", reference_level = "non_existent_group", verbose = FALSE),
        "reference_level 'non_existent_group' was not found in metadata\\$disease"
    )
})

test_that("run_diff_abundance issues warning when input abundance scale is not 'counts'", {
    d <- make_test_data(abundance_scale = "counts")
    validated <- suppressWarnings(
        validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
    )
    # Force scale label to 'percentage' to test MiCARA warning trigger
    validated$abundance_scale$final$taxa <- "percentage"

    expect_warning(
        withCallingHandlers(
            run_diff_abundance(validated, feature_type = "taxa", reference_level = "healthy", verbose = FALSE),
            warning = function(w) {
                if (grepl("< 3 categories|estimating sample-specific biases", w$message)) {
                    invokeRestart("muffleWarning")
                }
            }
        ),
        "ANCOM-BC2 log-linear models expect integer count data"
    )
})

test_that("run_diff_abundance warns and drops zero-variance fixed-effect covariates", {
    d <- make_test_data(abundance_scale = "counts")
    d$metadata$constant_cov <- 1 # zero variance column

    validated <- suppressWarnings(
        validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
    )

    expect_warning(
        withCallingHandlers(
            run_diff_abundance(validated, feature_type = "taxa", fix_effects = "constant_cov", reference_level = "healthy", verbose = FALSE),
            warning = function(w) {
                if (grepl("< 3 categories|estimating sample-specific biases", w$message)) {
                    invokeRestart("muffleWarning")
                }
            }
        ),
        "Covariate 'constant_cov' has zero variance"
    )
})

test_that("run_diff_abundance successfully fits model for single feature type ('taxa')", {
    skip_if_not_installed("ANCOMBC")
    skip_if_not_installed("phyloseq")

    d <- make_test_data(abundance_scale = "counts")
    validated <- suppressWarnings(
        validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
    )

    res <- suppressWarnings(
        run_diff_abundance(
            validated,
            feature_type = "taxa",
            reference_level = "healthy",
            verbose = FALSE
        )
    )

    expect_s3_class(res, "micara_diffab")
    expect_named(
        res,
        c("ancombc2_fit", "res", "significant", "feature_type", "fix_formula", "rand_formula", "reference_level", "use_robust"),
        ignore.order = TRUE
    )
    expect_equal(res$feature_type, "taxa")
    expect_true(is.data.frame(res$significant))
})

test_that("run_diff_abundance runs across both feature types ('both')", {
    skip_if_not_installed("ANCOMBC")
    skip_if_not_installed("phyloseq")

    d <- make_test_data(abundance_scale = "counts")
    validated <- suppressWarnings(
        validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
    )

    res <- suppressWarnings(
        run_diff_abundance(
            validated,
            feature_type = "both",
            reference_level = "healthy",
            verbose = FALSE
        )
    )

    expect_s3_class(res, "micara_diffab_list")
    expect_named(res, c("taxa", "pathways"))
    expect_s3_class(res$taxa, "micara_diffab")
    expect_s3_class(res$pathways, "micara_diffab")
})

test_that("run_diff_abundance sanitizes non-standard characters in disease labels", {
    skip_if_not_installed("ANCOMBC")
    skip_if_not_installed("phyloseq")

    d <- make_test_data(abundance_scale = "counts")
    d$metadata$disease <- ifelse(d$metadata$disease == "healthy", "Healthy Control", "T2D Stage-2")

    validated <- suppressWarnings(
        validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
    )

    res <- suppressWarnings(
        run_diff_abundance(
            validated,
            feature_type = "taxa",
            reference_level = "Healthy Control",
            verbose = FALSE
        )
    )

    expect_s3_class(res, "micara_diffab")
    expect_equal(res$reference_level, "Healthy.Control")
})
