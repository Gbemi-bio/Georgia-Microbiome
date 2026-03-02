#!/usr/bin/env Rscript
# =============================================================================
# Noise Assessment for ROOT Samples
# =============================================================================
# Rarefaction depth: 1000 reads per sample (CORRECTED)
# =============================================================================

# Load functions
source("noise_assessment_functions.R")

# Run analysis
results_root <- run_noise_assessment(
  phyloseq_file = "/Users/apple/Downloads/root_bac_sub_no_bad_Filtered_final_non_normalized (3).rds",
  sample_type = "root",
  rarefaction_depth = 1000,  # CORRECTED to match manuscript
  output_dir = ".",
  alpha_metrics = c("shannon", "simpson", "observed_species", "chao1", "pielou")
)

# Save results object
saveRDS(results_root, "root_noise_assessment_results.rds")
cat("\nResults object saved: root_noise_assessment_results.rds\n")
