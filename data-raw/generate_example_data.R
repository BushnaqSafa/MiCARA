
# ==============================================================================
# Script: generate_example_data.R
# Purpose: Generate synthetic example datasets with biological signal for MiCARA
# ==============================================================================

library(usethis)

# Set seed for reproducible synthetic data
set.seed(123) # Set seed for reproducible synthetic data

# 1. Metadata (30 samples)
example_metadata <- data.frame(
  sample_id    = paste0("Sample_", 1:30),
  disease      = rep(c("healthy", "Case"), each = 15),
  age_category = sample(c("Adult", "Senior"), 30, replace = TRUE),
  BMI          = round(runif(30, 18.5, 30.0), 1),
  gender       = sample(c("Male", "Female"), 30, replace = TRUE),
  study_name   = rep(c("Study_A", "Study_B"), times = 15),
  stringsAsFactors = FALSE
)
rownames(example_metadata) <- example_metadata$sample_id

# 2. Taxa Counts Matrix (20 taxa x 30 samples)
example_taxa <- as.data.frame(matrix(
  rpois(600, lambda = 50),
  nrow = 20,
  ncol = 30,
  dimnames = list(
    paste0("Taxon_", 1:20),
    example_metadata$sample_id
  )
))

# Inject signal for Taxon 1, 2, and 3 in Case group (Cols 16 to 30)
example_taxa["Taxon_1", 16:30] <- example_taxa["Taxon_1", 16:30] * 8
example_taxa["Taxon_2", 16:30] <- example_taxa["Taxon_2", 16:30] * 10
example_taxa["Taxon_3", 16:30] <- example_taxa["Taxon_3", 16:30] * 6

# 3. Pathway Counts Matrix (10 pathways x 30 samples)
example_pathways <- as.data.frame(matrix(
  rpois(300, lambda = 100),
  nrow = 10,
  ncol = 30,
  dimnames = list(
    paste0("Pathway_", 1:10),
    example_metadata$sample_id
  )
))

# Inject signal for Pathway 1 and 2 in Case group (Cols 16 to 30)
example_pathways["Pathway_1", 16:30] <- example_pathways["Pathway_1", 16:30] * 8
example_pathways["Pathway_2", 16:30] <- example_pathways["Pathway_2", 16:30] * 7

# Pathway_1 is made proportional to Taxon_1 plus a small random Poisson noise to create true correlation
example_pathways["Pathway_1", ] <- example_taxa["Taxon_1", ] * 2 + rpois(30, lambda = 5)


# Save binary .rda objects to the package's data/ folder
usethis::use_data(example_taxa, overwrite = TRUE)
usethis::use_data(example_pathways, overwrite = TRUE)
usethis::use_data(example_metadata, overwrite = TRUE)

message("✅ Example datasets successfully generated with injected differential signal!")

