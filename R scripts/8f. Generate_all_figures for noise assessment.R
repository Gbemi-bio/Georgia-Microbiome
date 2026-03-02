#!/usr/bin/env Rscript
# =============================================================================
# Generate All Figures for Noise Assessment
# =============================================================================
# Creates publication-quality figures for all sample types
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Set theme
theme_set(theme_bw(base_size = 12))

# Define sample types and their colors
sample_types <- c("root", "soil", "stalk", "rhizosphere")
sample_colors <- c(
  "root" = "#E78AC3",
  "soil" = "#8DA0CB",
  "stalk" = "#66C2A5",
  "rhizosphere" = "#FC8D62"
)

cat("=============================================================================\n")
cat("GENERATING ALL FIGURES\n")
cat("=============================================================================\n\n")

# =============================================================================
# Function: Create Alpha Diversity CV Plot
# =============================================================================

create_alpha_cv_plot <- function(sample_type) {
  alpha_summary <- read.csv(paste0(sample_type, "_alpha_noise_summary.csv"))

  cv_data <- alpha_summary %>%
    select(Metric, Mean_CV, SD_CV) %>%
    mutate(
      Metric = factor(Metric,
                      levels = c("shannon", "simpson", "pielou", "observed_species", "chao1"),
                      labels = c("Shannon", "Simpson", "Pielou", "Observed\nSpecies", "Chao1"))
    )

  p <- ggplot(cv_data, aes(x = Metric, y = Mean_CV)) +
    geom_bar(stat = "identity", fill = sample_colors[sample_type],
             color = "black", alpha = 0.8) +
    geom_errorbar(aes(ymin = Mean_CV, ymax = Mean_CV + SD_CV),
                  width = 0.3, linewidth = 0.5) +
    labs(
      title = paste("Alpha Diversity CV -", toupper(sample_type)),
      x = "Diversity Metric",
      y = "Coefficient of Variation (%)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    ) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red", alpha = 0.5)

  return(p)
}

# =============================================================================
# Function: Create Beta Diversity Comparison Plot
# =============================================================================

create_beta_comparison_plot <- function(sample_type) {
  beta_summary <- read.csv(paste0(sample_type, "_beta_noise_summary.csv"))
  paired_dist <- read.csv(paste0(sample_type, "_beta_noise_paired_distances.csv"))

  # Create data for plotting
  beta_data <- paired_dist %>%
    select(FieldID, Growing_Degree_Days, within_pair_distance) %>%
    mutate(Type = "Within-Pair")

  # Generate random sample
  set.seed(123)
  n_random <- nrow(beta_data)
  random_mean <- beta_summary$Mean[beta_summary$Type == "Random_Pair"]
  random_sd <- beta_summary$SD[beta_summary$Type == "Random_Pair"]

  # Handle NA cases (e.g., when only 1 pair)
  if (!is.na(random_sd)) {
    random_data <- data.frame(
      FieldID = NA,
      Growing_Degree_Days = NA,
      within_pair_distance = rnorm(n_random, mean = random_mean, sd = random_sd),
      Type = "Random Pairs"
    )
  } else {
    random_data <- data.frame(
      FieldID = NA,
      Growing_Degree_Days = NA,
      within_pair_distance = random_mean,
      Type = "Random Pairs"
    )
  }

  combined_beta <- bind_rows(beta_data, random_data) %>%
    mutate(Type = factor(Type, levels = c("Within-Pair", "Random Pairs")))

  p <- ggplot(combined_beta, aes(x = Type, y = within_pair_distance, fill = Type)) +
    geom_violin(alpha = 0.6, color = "black") +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(data = subset(combined_beta, Type == "Within-Pair"),
                width = 0.1, alpha = 0.4, size = 1.5) +
    scale_fill_manual(values = c("Within-Pair" = sample_colors[sample_type],
                                  "Random Pairs" = "gray70")) +
    labs(
      title = paste("Beta Diversity -", toupper(sample_type)),
      x = "Comparison Type",
      y = "Weighted UniFrac Distance"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  return(p)
}

# =============================================================================
# Function: Create Paired Connection Plot (NEW - Aduragbemi's request)
# =============================================================================

create_paired_connection_plot <- function(sample_type) {
  # Check if PCoA files exist
  pcoa_file <- paste0(sample_type, "_pcoa_coordinates.csv")
  conn_file <- paste0(sample_type, "_pcoa_connections.csv")

  if (!file.exists(pcoa_file) || !file.exists(conn_file)) {
    cat("Skipping paired connection plot for", sample_type, "(insufficient samples)\n")
    return(NULL)
  }

  pcoa_df <- read.csv(pcoa_file)
  connections <- read.csv(conn_file)

  # Load results to get variance explained
  results <- readRDS(paste0(sample_type, "_noise_assessment_results.rds"))
  var_exp <- results$variance_explained

  p <- ggplot() +
    # Draw connection lines first (so they're behind points)
    geom_segment(data = connections,
                 aes(x = PC1_p1, y = PC2_p1, xend = PC1_p2, yend = PC2_p2),
                 color = "gray50", alpha = 0.5, linewidth = 0.5) +
    # Draw points (p1 and p2 use the same shape since assignment is arbitrary)
    geom_point(data = pcoa_df,
               aes(x = PC1, y = PC2, color = factor(Growing_Degree_Days)),
               shape = 16, size = 8, alpha = 0.8) +
    scale_color_manual(values = c("600" = "#E78AC3", "1400" = "#8DA0CB"),
                       name = "GDD") +
    labs(
      title = paste("Paired Sample Clustering -", toupper(sample_type)),
      subtitle = "Lines connect paired samples from the same location",
      x = paste0("PC1 (", var_exp[1], "%)"),
      y = paste0("PC2 (", var_exp[2], "%)")
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "right"
    )

  return(p)
}

# =============================================================================
# Generate Figures for Each Sample Type
# =============================================================================

for (sample_type in sample_types) {
  cat("\n--- Generating figures for", toupper(sample_type), "---\n")

  # Alpha CV plot
  p_alpha <- create_alpha_cv_plot(sample_type)
  ggsave(paste0("Figure_", sample_type, "_alpha_CV.png"), p_alpha,
         width = 8, height = 6, dpi = 300)
  ggsave(paste0("Figure_", sample_type, "_alpha_CV.pdf"), p_alpha,
         width = 8, height = 6)
  cat("  Created: Figure_", sample_type, "_alpha_CV.png/pdf\n", sep = "")

  # Beta comparison plot
  p_beta <- create_beta_comparison_plot(sample_type)
  ggsave(paste0("Figure_", sample_type, "_beta_comparison.png"), p_beta,
         width = 8, height = 6, dpi = 300)
  ggsave(paste0("Figure_", sample_type, "_beta_comparison.pdf"), p_beta,
         width = 8, height = 6)
  cat("  Created: Figure_", sample_type, "_beta_comparison.png/pdf\n", sep = "")

  # Paired connection plot (NEW)
  p_paired <- create_paired_connection_plot(sample_type)
  if (!is.null(p_paired)) {
    ggsave(paste0("Figure_", sample_type, "_paired_connections.png"), p_paired,
           width = 10, height = 8, dpi = 300)
    ggsave(paste0("Figure_", sample_type, "_paired_connections.pdf"), p_paired,
           width = 10, height = 8)
    cat("  Created: Figure_", sample_type, "_paired_connections.png/pdf\n", sep = "")
  }
}

# =============================================================================
# Create Combined Comparison Figure (All Sample Types)
# =============================================================================

cat("\n--- Creating combined comparison figures ---\n")

# Combined Alpha CV
alpha_combined <- bind_rows(lapply(sample_types, function(st) {
  read.csv(paste0(st, "_alpha_noise_summary.csv")) %>%
    mutate(Sample_Type = st)
}))

alpha_combined <- alpha_combined %>%
  mutate(
    Metric = factor(Metric,
                    levels = c("shannon", "simpson", "pielou", "observed_species", "chao1"),
                    labels = c("Shannon", "Simpson", "Pielou", "Observed Species", "Chao1")),
    Sample_Type = factor(Sample_Type,
                         levels = sample_types,
                         labels = c("Root", "Soil", "Stalk", "Rhizosphere"))
  )

p_combined_alpha <- ggplot(alpha_combined, aes(x = Metric, y = Mean_CV, fill = Sample_Type)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", alpha = 0.8) +
  geom_errorbar(aes(ymin = Mean_CV, ymax = Mean_CV + SD_CV),
                position = position_dodge(0.9), width = 0.25) +
  scale_fill_manual(values = sample_colors) +
  labs(
    title = "Alpha Diversity Noise Across All Sample Types",
    x = "Diversity Metric",
    y = "Mean Coefficient of Variation (%)",
    fill = "Sample Type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red", alpha = 0.5)

ggsave("Figure_Combined_Alpha_CV.png", p_combined_alpha,
       width = 12, height = 8, dpi = 300)
ggsave("Figure_Combined_Alpha_CV.pdf", p_combined_alpha,
       width = 12, height = 8)
cat("Created: Figure_Combined_Alpha_CV.png/pdf\n")

# Combined Beta Diversity
beta_combined <- bind_rows(lapply(sample_types, function(st) {
  read.csv(paste0(st, "_beta_noise_summary.csv")) %>%
    filter(Type == "Within_Pair") %>%
    mutate(Sample_Type = st)
}))

beta_combined <- beta_combined %>%
  mutate(
    Sample_Type = factor(Sample_Type,
                         levels = sample_types,
                         labels = c("Root", "Soil", "Stalk", "Rhizosphere"))
  )

p_combined_beta <- ggplot(beta_combined, aes(x = Sample_Type, y = Mean, fill = Sample_Type)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  geom_errorbar(aes(ymin = Mean, ymax = Mean + SD),
                width = 0.3, linewidth = 0.5) +
  scale_fill_manual(values = sample_colors) +
  labs(
    title = "Within-Pair Beta Diversity Across Sample Types",
    x = "Sample Type",
    y = "Mean Weighted UniFrac Distance"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "none"
  )

ggsave("Figure_Combined_Beta_Distance.png", p_combined_beta,
       width = 8, height = 6, dpi = 300)
ggsave("Figure_Combined_Beta_Distance.pdf", p_combined_beta,
       width = 8, height = 6)
cat("Created: Figure_Combined_Beta_Distance.png/pdf\n")

cat("\n=============================================================================\n")
cat("ALL FIGURES GENERATED SUCCESSFULLY\n")
cat("=============================================================================\n\n")

cat("Individual figures created for each sample type:\n")
cat("  - Alpha CV plots\n")
cat("  - Beta diversity comparison plots\n")
cat("  - Paired connection plots (PCoA with connecting lines)\n\n")
cat("Combined figures:\n")
cat("  - Combined alpha CV comparison\n")
cat("  - Combined beta diversity comparison\n\n")

cat("Total figures created:\n")
cat("  PNG files:", length(list.files(pattern = "Figure_.*\\.png$")), "\n")
cat("  PDF files:", length(list.files(pattern = "Figure_.*\\.pdf$")), "\n")
