# Helper to create lightweight mock micara_interactions objects for testing
make_mock_sankey_interactions <- function() {
  # DiseaseA has mixed positive and negative links
  links_a <- data.frame(
    taxon = c("Taxon_1", "Taxon_2", "Taxon_3", "Taxon_4"),
    pathway = c("Pathway_1", "Pathway_1", "Pathway_2", "Pathway_2"),
    rho = c(0.45, -0.35, 0.60, -0.25),
    p = c(0.001, 0.01, 0.0005, 0.02),
    q = c(0.005, 0.02, 0.002, 0.04),
    direction = c("positive", "negative", "positive", "negative"),
    stringsAsFactors = FALSE
  )
  
  node_directions_a <- data.frame(
    name = c("Taxon_1", "Taxon_2", "Taxon_3", "Taxon_4", "Pathway_1", "Pathway_2"),
    type = c(rep("Taxon", 4), rep("Pathway", 2)),
    group = c("up", "down", "up", "down", "up", "unchanged"),
    stringsAsFactors = FALSE
  )
  
  disease_a <- list(
    status = "analysed",
    links = links_a,
    node_directions = node_directions_a
  )
  
  # DiseaseB has ONLY positive links
  links_b <- data.frame(
    taxon = c("Taxon_1", "Taxon_2"),
    pathway = c("Pathway_1", "Pathway_1"),
    rho = c(0.50, 0.40),
    p = c(0.001, 0.005),
    q = c(0.01, 0.02),
    direction = c("positive", "positive"),
    stringsAsFactors = FALSE
  )
  
  node_directions_b <- data.frame(
    name = c("Taxon_1", "Taxon_2", "Pathway_1"),
    type = c("Taxon", "Taxon", "Pathway"),
    group = c("up", "up", "down"),
    stringsAsFactors = FALSE
  )
  
  disease_b <- list(
    status = "analysed",
    links = links_b,
    node_directions = node_directions_b
  )
  
  # DiseaseEmpty has 0 links passing thresholds
  disease_empty <- list(
    status = "skipped",
    links = links_a[0, ],
    node_directions = node_directions_a[0, ]
  )
  
  res <- list(
    DiseaseA = disease_a,
    Disease.Sanitized = disease_b,
    DiseaseEmpty = disease_empty
  )
  
  class(res) <- c("micara_interactions", "list")
  res
}

test_that(".to_js_array creates valid JavaScript array strings", {
  js_str <- MiCARA:::.to_js_array(c("up", "down", "unchanged"))
  expect_equal(js_str, "[\"up\", \"down\", \"unchanged\"]")
})

test_that("plot_interaction_sankey validates input object class and parameters", {
  skip_if_not_installed("networkD3")
  interactions <- make_mock_sankey_interactions()
  
  # Invalid S3 class
  expect_error(
    plot_interaction_sankey(list(), disease = "DiseaseA"),
    "'interactions' must be the output of compute_residualCLR_correlations"
  )
  
  # Invalid disease inputs (non-scalar / empty / NA)
  expect_error(
    plot_interaction_sankey(interactions, disease = c("DiseaseA", "DiseaseB")),
    "'disease' must be a single non-empty character string"
  )
  expect_error(
    plot_interaction_sankey(interactions, disease = NA_character_),
    "'disease' must be a single non-empty character string"
  )
  
  # Invalid direction parameter (match.arg handling)
  expect_error(
    plot_interaction_sankey(interactions, disease = "DiseaseA", direction = "invalid_dir")
  )
  
  # Disease not present in interactions object
  expect_error(
    plot_interaction_sankey(interactions, disease = "NonExistentDisease"),
    "'NonExistentDisease' not found in interactions. Available diseases: DiseaseA, Disease.Sanitized, DiseaseEmpty"
  )
  
  # Unnamed node_colors vector
  expect_error(
    plot_interaction_sankey(interactions, disease = "DiseaseA", node_colors = c("#2196F3", "#F44336")),
    "'node_colors' must be a named vector"
  )
  
  # Invalid top_n_per_pathway parameter
  expect_error(
    plot_interaction_sankey(interactions, disease = "DiseaseA", top_n_per_pathway = 0),
    "'top_n_per_pathway' must be a positive integer"
  )
})

test_that("plot_interaction_sankey handles empty link tables and missing directions", {
  skip_if_not_installed("networkD3")
  interactions <- make_mock_sankey_interactions()
  
  # Disease has zero total links
  expect_error(
    plot_interaction_sankey(interactions, disease = "DiseaseEmpty"),
    "DiseaseEmpty has no links passing the significance/strength cutoffs; nothing to plot."
  )
  
  # Disease has no negative links when requested
  expect_error(
    plot_interaction_sankey(interactions, disease = "Disease.Sanitized", direction = "negative"),
    "Disease.Sanitized has no negative links passing the cutoffs; nothing to plot."
  )
})

test_that("plot_interaction_sankey resolves non-syntactic or sanitized disease names", {
  skip_if_not_installed("networkD3")
  interactions <- make_mock_sankey_interactions()
  
  # Requesting "Disease-Sanitized" resolves to "Disease.Sanitized"
  widget <- plot_interaction_sankey(interactions, disease = "Disease-Sanitized")
  expect_s3_class(widget, "sankeyNetwork")
  expect_s3_class(widget, "htmlwidget")
})

test_that("plot_interaction_sankey filters by direction accurately", {
  skip_if_not_installed("networkD3")
  interactions <- make_mock_sankey_interactions()
  
  # Both directions (default)
  widget_both <- plot_interaction_sankey(interactions, disease = "DiseaseA", direction = "both")
  expect_equal(nrow(widget_both$x$links), 4)
  
  # Positive direction only
  widget_pos <- plot_interaction_sankey(interactions, disease = "DiseaseA", direction = "positive")
  expect_equal(nrow(widget_pos$x$links), 2)
  
  # Negative direction only
  widget_neg <- plot_interaction_sankey(interactions, disease = "DiseaseA", direction = "negative")
  expect_equal(nrow(widget_neg$x$links), 2)
})

test_that("plot_interaction_sankey trims links using top_n_per_pathway", {
  skip_if_not_installed("networkD3")
  interactions <- make_mock_sankey_interactions()
  
  # Top 1 per pathway should select only the single highest abs(rho) per pathway (2 total)
  widget_top1 <- plot_interaction_sankey(
    interactions,
    disease = "DiseaseA",
    top_n_per_pathway = 1
  )
  
  expect_equal(nrow(widget_top1$x$links), 2)
})

test_that("plot_interaction_sankey renders networkD3 widget and handles missing node_directions", {
  skip_if_not_installed("networkD3")
  interactions <- make_mock_sankey_interactions()
  
  custom_colors <- c(up = "blue", down = "red", unchanged = "grey", stable = "grey")
  widget <- plot_interaction_sankey(
    interactions,
    disease = "DiseaseA",
    node_colors = custom_colors,
    font_size = 14,
    node_width = 40
  )
  
  expect_s3_class(widget, "sankeyNetwork")
  expect_equal(widget$x$options$fontSize, 14)
  expect_equal(widget$x$options$nodeWidth, 40)
  
  # Test fallback when node_directions is missing
  interactions$DiseaseA$node_directions <- NULL
  widget_fallback <- plot_interaction_sankey(interactions, disease = "DiseaseA")
  expect_true(all(widget_fallback$x$nodes$group == "stable"))
})

test_that("plot_interaction_sankey exports files and strips duplicate extension cleanly", {
  skip_if_not_installed("networkD3")
  skip_if_not_installed("htmlwidgets")
  interactions <- make_mock_sankey_interactions()
  
  temp_prefix <- tempfile("sankey_test_")
  temp_with_ext <- paste0(temp_prefix, ".html")
  
  # Pass path WITH .html extension to test sub() extension stripping
  res <- suppressWarnings(
    plot_interaction_sankey(
      interactions,
      disease = "DiseaseA",
      save_path = temp_with_ext
    )
  )
  
  expect_s3_class(res, "sankeyNetwork")
  expect_true(file.exists(temp_with_ext))
  expect_false(file.exists(paste0(temp_with_ext, ".html")))
  
  # Cleanup
  if (file.exists(temp_with_ext)) unlink(temp_with_ext)
  png_expected <- paste0(temp_prefix, ".png")
  if (file.exists(png_expected)) unlink(png_expected)
})