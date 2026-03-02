#!/usr/bin/env Rscript
# =============================================================================
# Noise Assessment for RHIZOSPHERE Samples
# =============================================================================
# Rarefaction depth: 3000 reads per sample
# =============================================================================

# Load functions
source("noise_assessment_functions.R")

# Run analysis
results_rhizosphere <- run_noise_assessment(
  phyloseq_file = "/Users/apple/Downloads/rhizo_bac_sub_no_bad_Filtered_final_non_normalized (1).rds",
  sample_type = "rhizosphere",
  rarefaction_depth = 3000,
  output_dir = ".",
  alpha_metrics = c("shannon", "simpson", "observed_species", "chao1", "pielou")
)

# Save results object
saveRDS(results_rhizosphere, "rhizosphere_noise_assessment_results.rds")
cat("\nResults object saved: rhizosphere_noise_assessment_results.rds\n")
