#!/usr/bin/env Rscript
# =============================================================================
# Noise Assessment Functions - Universal Framework
# =============================================================================
# Reusable functions for noise assessment across all sample types
# =============================================================================

library(phyloseq)
library(ape)
library(MicrobiomeStat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)

#' Universal Noise Assessment Function
#'
#' @param phyloseq_file Path to phyloseq RDS file
#' @param sample_type Name of sample type (e.g., "root", "soil")
#' @param rarefaction_depth Rarefaction depth for this sample type
#' @param output_dir Directory to save results
#' @param alpha_metrics Vector of alpha diversity metrics to calculate
#'
#' @return List containing all results and file paths
run_noise_assessment <- function(phyloseq_file,
                                  sample_type,
                                  rarefaction_depth,
                                  output_dir = ".",
                                  alpha_metrics = c("shannon", "simpson", "observed_species", "chao1", "pielou")) {

  cat("\n")
  cat("=============================================================================\n")
  cat("NOISE ASSESSMENT FOR", toupper(sample_type), "SAMPLES\n")
  cat("Rarefaction depth:", rarefaction_depth, "reads\n")
  cat("=============================================================================\n\n")

  # -------------------------------------------------------------------------
  # 1. Load and Prepare Data
  # -------------------------------------------------------------------------

  cat("Loading data from:", phyloseq_file, "\n")
  ps_data <- readRDS(phyloseq_file)

  cat("Original samples:", nsamples(ps_data), "\n")
  cat("Original taxa:", ntaxa(ps_data), "\n\n")

  # Check and fix tree
  cat("Checking phylogenetic tree...\n")
  if (!is.null(phy_tree(ps_data, errorIfNULL = FALSE))) {
    ps_tree <- phy_tree(ps_data)
    if(!is.binary(ps_tree)) {
      cat("Converting multifurcating tree to binary...\n")
      phy_tree(ps_data) <- multi2di(ps_tree)
    }
  } else {
    cat("Warning: No phylogenetic tree found. Faith PD will not be calculated.\n")
  }

  # -------------------------------------------------------------------------
  # 2. Filter to Complete Pairs
  # -------------------------------------------------------------------------

  cat("\nFiltering to complete pairs...\n")
  meta <- data.frame(sample_data(ps_data))

  # Identify complete pairs
  complete_pairs <- meta %>%
    group_by(FieldID, Growing_Degree_Days) %>%
    summarise(
      has_both = sum(time == "p1") >= 1 & sum(time == "p2") >= 1,
      .groups = "drop"
    ) %>%
    filter(has_both)

  cat("Complete pairs found:", nrow(complete_pairs), "\n\n")

  # Filter samples
  samples_to_keep <- meta %>%
    inner_join(complete_pairs, by = c("FieldID", "Growing_Degree_Days")) %>%
    pull(Samples)

  ps_paired <- prune_samples(samples_to_keep, ps_data)
  cat("Samples after filtering to complete pairs:", nsamples(ps_paired), "\n")

  # -------------------------------------------------------------------------
  # 3. Rarefy Data
  # -------------------------------------------------------------------------

  cat("\nRarefying to", rarefaction_depth, "reads per sample...\n")
  set.seed(123)
  ps_rarefied <- rarefy_even_depth(ps_paired,
                                    sample.size = rarefaction_depth,
                                    replace = FALSE)
  cat("Samples after rarefaction:", nsamples(ps_rarefied), "\n\n")

  # -------------------------------------------------------------------------
  # 4. Alpha Diversity Noise Assessment
  # -------------------------------------------------------------------------

  cat("=============================================================================\n")
  cat("ALPHA DIVERSITY NOISE ASSESSMENT\n")
  cat("=============================================================================\n\n")

  # Convert to MicrobiomeStat object
  data.obj <- mStat_convert_phyloseq_to_data_obj(ps_rarefied)

  # Calculate alpha diversity
  cat("Calculating alpha diversity metrics...\n")
  alpha_data <- mStat_calculate_alpha_diversity(
    x = data.obj$feature.tab,
    alpha.name = alpha_metrics
  )

  # Combine with metadata
  alpha_df <- bind_cols(alpha_data) %>%
    rownames_to_column("sample") %>%
    inner_join(
      data.obj$meta.dat %>% rownames_to_column("sample"),
      by = "sample"
    )

  # Calculate noise metrics
  cat("\nCalculating alpha diversity noise for each location-GDD pair...\n")

  alpha_noise_results <- alpha_df %>%
    group_by(FieldID, Growing_Degree_Days) %>%
    filter(n() == 2, all(c("p1", "p2") %in% time)) %>%
    summarise(
      across(
        all_of(alpha_metrics),
        list(
          p1 = ~.[time == "p1"][1],
          p2 = ~.[time == "p2"][1],
          abs_diff = ~abs(.[time == "p2"][1] - .[time == "p1"][1]),
          rel_diff_pct = ~abs((.[time == "p2"][1] - .[time == "p1"][1]) / .[time == "p1"][1] * 100),
          cv = ~sd(c(.[time == "p1"][1], .[time == "p2"][1])) / mean(c(.[time == "p1"][1], .[time == "p2"][1])) * 100
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )

  # Summary statistics
  cat("\n--- Summary of Alpha Diversity Noise Across All Pairs ---\n\n")

  alpha_noise_summary <- data.frame(
    Metric = alpha_metrics,
    Mean_Abs_Diff = sapply(alpha_metrics, function(m) mean(alpha_noise_results[[paste0(m, "_abs_diff")]], na.rm = TRUE)),
    SD_Abs_Diff = sapply(alpha_metrics, function(m) sd(alpha_noise_results[[paste0(m, "_abs_diff")]], na.rm = TRUE)),
    Mean_Rel_Diff_Pct = sapply(alpha_metrics, function(m) mean(alpha_noise_results[[paste0(m, "_rel_diff_pct")]], na.rm = TRUE)),
    SD_Rel_Diff_Pct = sapply(alpha_metrics, function(m) sd(alpha_noise_results[[paste0(m, "_rel_diff_pct")]], na.rm = TRUE)),
    Mean_CV = sapply(alpha_metrics, function(m) mean(alpha_noise_results[[paste0(m, "_cv")]], na.rm = TRUE)),
    SD_CV = sapply(alpha_metrics, function(m) sd(alpha_noise_results[[paste0(m, "_cv")]], na.rm = TRUE))
  )

  print(alpha_noise_summary)

  # -------------------------------------------------------------------------
  # 5. Beta Diversity Noise Assessment
  # -------------------------------------------------------------------------

  cat("\n\n=============================================================================\n")
  cat("BETA DIVERSITY NOISE ASSESSMENT\n")
  cat("=============================================================================\n\n")

  # Calculate weighted UniFrac distances
  cat("Calculating Weighted UniFrac distances...\n")
  wunifrac_dist <- phyloseq::distance(ps_rarefied, method = "wunifrac")
  wunifrac_mat <- as.matrix(wunifrac_dist)

  # Extract paired distances
  cat("Extracting within-pair distances...\n")

  meta_rarefied <- data.frame(sample_data(ps_rarefied))
  meta_rarefied$sample_name <- rownames(meta_rarefied)

  paired_distances <- meta_rarefied %>%
    group_by(FieldID, Growing_Degree_Days) %>%
    filter(n() == 2, all(c("p1", "p2") %in% time)) %>%
    summarise(
      sample_p1 = sample_name[time == "p1"][1],
      sample_p2 = sample_name[time == "p2"][1],
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      within_pair_distance = wunifrac_mat[sample_p1, sample_p2]
    ) %>%
    ungroup()

  cat("\n--- Summary of Beta Diversity Noise (Weighted UniFrac) ---\n\n")
  cat("Number of complete pairs:", nrow(paired_distances), "\n")
  cat("Mean within-pair distance:", round(mean(paired_distances$within_pair_distance), 4), "\n")
  cat("SD within-pair distance:", round(sd(paired_distances$within_pair_distance), 4), "\n")
  cat("Median within-pair distance:", round(median(paired_distances$within_pair_distance), 4), "\n")

  # Compare with random pairs
  set.seed(123)
  n_random <- min(1000, length(wunifrac_dist))
  random_pairs <- sample(wunifrac_dist, n_random, replace = FALSE)

  cat("\n--- Comparison with Random Pairs ---\n\n")
  cat("Mean random-pair distance:", round(mean(random_pairs), 4), "\n")
  cat("SD random-pair distance:", round(sd(random_pairs), 4), "\n")

  # Statistical test
  wilcox_test <- wilcox.test(paired_distances$within_pair_distance, random_pairs, alternative = "less")
  cat("\nWilcoxon test (paired < random):\n")
  cat("  W =", wilcox_test$statistic, "\n")
  cat("  p-value =", format(wilcox_test$p.value, scientific = TRUE), "\n\n")

  # -------------------------------------------------------------------------
  # 6. PCoA for Paired Connection Plot
  # -------------------------------------------------------------------------

  cat("=============================================================================\n")
  cat("CALCULATING PCOA FOR PAIRED CONNECTION PLOT\n")
  cat("=============================================================================\n\n")

  # Check if we have enough samples for PCoA
  n_samples <- nsamples(ps_rarefied)

  if (n_samples >= 3) {
    # Perform PCoA
    pcoa_result <- cmdscale(wunifrac_mat, k = 2, eig = TRUE)

    # Calculate variance explained
    variance_explained <- round(100 * pcoa_result$eig[1:2] / sum(abs(pcoa_result$eig)), 2)

    # Create PCoA dataframe
    pcoa_df <- data.frame(
      sample = rownames(pcoa_result$points),
      PC1 = pcoa_result$points[, 1],
      PC2 = pcoa_result$points[, 2]
    ) %>%
      left_join(meta_rarefied, by = c("sample" = "sample_name"))

    # Create connection data
    connections <- paired_distances %>%
      left_join(pcoa_df %>% select(sample, PC1_p1 = PC1, PC2_p1 = PC2),
                by = c("sample_p1" = "sample")) %>%
      left_join(pcoa_df %>% select(sample, PC1_p2 = PC1, PC2_p2 = PC2),
                by = c("sample_p2" = "sample"))

    cat("PCoA calculated. Variance explained:\n")
    cat("  PC1:", variance_explained[1], "%\n")
    cat("  PC2:", variance_explained[2], "%\n\n")
  } else {
    cat("WARNING: Not enough samples (n =", n_samples, ") for PCoA calculation.\n")
    cat("Need at least 3 samples. PCoA will be skipped.\n\n")
    pcoa_df <- NULL
    connections <- NULL
    variance_explained <- NULL
  }

  # -------------------------------------------------------------------------
  # 7. Save Results
  # -------------------------------------------------------------------------

  cat("=============================================================================\n")
  cat("SAVING RESULTS\n")
  cat("=============================================================================\n\n")

  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Define file prefix
  prefix <- file.path(output_dir, paste0(sample_type, "_"))

  # Save alpha diversity noise
  write.csv(alpha_noise_results,
            paste0(prefix, "alpha_noise_by_location.csv"),
            row.names = FALSE)
  cat("Saved:", paste0(prefix, "alpha_noise_by_location.csv\n"))

  write.csv(alpha_noise_summary,
            paste0(prefix, "alpha_noise_summary.csv"),
            row.names = FALSE)
  cat("Saved:", paste0(prefix, "alpha_noise_summary.csv\n"))

  # Save beta diversity noise
  write.csv(paired_distances,
            paste0(prefix, "beta_noise_paired_distances.csv"),
            row.names = FALSE)
  cat("Saved:", paste0(prefix, "beta_noise_paired_distances.csv\n"))

  # Beta summary
  beta_summary <- data.frame(
    Metric = "Weighted_UniFrac",
    Type = c("Within_Pair", "Random_Pair"),
    Mean = c(mean(paired_distances$within_pair_distance), mean(random_pairs)),
    SD = c(sd(paired_distances$within_pair_distance), sd(random_pairs)),
    Median = c(median(paired_distances$within_pair_distance), median(random_pairs)),
    Min = c(min(paired_distances$within_pair_distance), min(random_pairs)),
    Max = c(max(paired_distances$within_pair_distance), max(random_pairs)),
    Q25 = c(quantile(paired_distances$within_pair_distance, 0.25), quantile(random_pairs, 0.25)),
    Q75 = c(quantile(paired_distances$within_pair_distance, 0.75), quantile(random_pairs, 0.75))
  )

  write.csv(beta_summary,
            paste0(prefix, "beta_noise_summary.csv"),
            row.names = FALSE)
  cat("Saved:", paste0(prefix, "beta_noise_summary.csv\n"))

  # Save PCoA data (only if calculated)
  if (!is.null(pcoa_df)) {
    write.csv(pcoa_df,
              paste0(prefix, "pcoa_coordinates.csv"),
              row.names = FALSE)
    cat("Saved:", paste0(prefix, "pcoa_coordinates.csv\n"))

    write.csv(connections,
              paste0(prefix, "pcoa_connections.csv"),
              row.names = FALSE)
    cat("Saved:", paste0(prefix, "pcoa_connections.csv\n"))
  } else {
    cat("PCoA data not saved (insufficient samples)\n")
  }

  cat("\n=============================================================================\n")
  cat("ANALYSIS COMPLETE FOR", toupper(sample_type), "\n")
  cat("=============================================================================\n\n")

  # Return results
  return(list(
    sample_type = sample_type,
    rarefaction_depth = rarefaction_depth,
    n_pairs = nrow(paired_distances),
    alpha_summary = alpha_noise_summary,
    alpha_by_location = alpha_noise_results,
    beta_summary = beta_summary,
    paired_distances = paired_distances,
    pcoa_df = pcoa_df,
    connections = connections,
    variance_explained = variance_explained,
    wilcox_test = wilcox_test,
    random_pairs = random_pairs
  ))
}

cat("Noise assessment functions loaded successfully.\n")
cat("Use: results <- run_noise_assessment(file, sample_type, rarefaction_depth)\n")
