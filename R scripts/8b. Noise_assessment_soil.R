#!/usr/bin/env Rscript
# =============================================================================
# Noise Assessment for SOIL Samples
# =============================================================================
# Rarefaction depth: 2000 reads per sample
# =============================================================================

# Load functions
source("noise_assessment_functions.R")

# Run analysis
results_soil <- run_noise_assessment(
  phyloseq_file = "/Users/apple/Downloads/soil_bac_sub_no_bad_Filtered_final_non_normalized.rds",
  sample_type = "soil",
  rarefaction_depth = 2000,
  output_dir = ".",
  alpha_metrics = c("shannon", "simpson", "observed_species", "chao1", "pielou")
)

# Save results object
saveRDS(results_soil, "soil_noise_assessment_results.rds")
cat("\nResults object saved: soil_noise_assessment_results.rds\n")
