#!/usr/bin/env Rscript
# =============================================================================
# Noise Assessment for STALK Samples
# =============================================================================
# Rarefaction depth: 500 reads per sample
# =============================================================================

# Load functions
source("noise_assessment_functions.R")

# Run analysis
results_stalk <- run_noise_assessment(
  phyloseq_file = "/Users/apple/Downloads/stalk_bac_sub_no_bad_Filtered_final_non_normalized.rds",
  sample_type = "stalk",
  rarefaction_depth = 500,
  output_dir = ".",
  alpha_metrics = c("shannon", "simpson", "observed_species", "chao1", "pielou")
)

# Save results object
saveRDS(results_stalk, "stalk_noise_assessment_results.rds")
cat("\nResults object saved: stalk_noise_assessment_results.rds\n")
