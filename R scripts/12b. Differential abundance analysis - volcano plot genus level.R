# Load required libraries
library(ggplot2)
library(ggrepel)

#read in the data
res_sig_soil_1 <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/res_sig_soil_1.rds")
res_sig_rhizo_1 <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/res_sig_rhizo_1.rds")
res_sig_root_1 <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/res_sig_root_1.rds")
res_sig_stalk_1 <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/res_sig_stalk_1.rds")


# ========================================================================
# Volcano plot of Differential abundance analysis of sample types by Genus 
# ========================================================================

# For Soil sample

# Create -log10 pvalue safely
res_sig_soil_1$logP <- -log10(res_sig_soil_1$pvalue)

# Remove infinite values if any
res_sig_soil_1 <- res_sig_soil_1[is.finite(res_sig_soil_1$logP), ]

# Define differential expression category
res_sig_soil_1$diffexpressed <- "Not Significant"

res_sig_soil_1$diffexpressed[
  res_sig_soil_1$log2FoldChange > 0.6 & res_sig_soil_1$pvalue < 1e-6
] <- "GDD_1400"

res_sig_soil_1$diffexpressed[
  res_sig_soil_1$log2FoldChange < -0.6 & res_sig_soil_1$pvalue < 1e-6
] <- "GDD_600"

# Label ONLY significant taxa
res_sig_soil_1$delabel <- ifelse(
  res_sig_soil_1$diffexpressed != "Not Significant",
  res_sig_soil_1$Genus,
  NA
)

# Volcano plot
soil_volcano <- ggplot(res_sig_soil_1,
                       aes(x = log2FoldChange,
                           y = logP,
                           color = diffexpressed)) +
  geom_point(size = 15, alpha = 0.9) +
  
  geom_vline(xintercept = c(-0.6, 0.6),
             color = "#2C6BED",
             linewidth = 0.8) +
  
  geom_hline(yintercept = -log10(1e-6),
             color = "#2C6BED",
             linewidth = 0.8) +
  
  geom_text_repel(
    data = subset(res_sig_soil_1, !is.na(delabel)),
    aes(label = delabel),
    size = 4.5,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  
  scale_color_manual(
    values = c(
      "GDD_600" = "#2C6BED",      
      "Not Significant" = "grey70",
      "GDD_1400" = "#D62728"      
    ),
    name = "More Abundant In:"
  ) +
  
  labs(
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~italic(p)~value)
  ) +
  
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text  = element_text(size = 14),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank()
  )

soil_volcano


#save plot
ggsave("soil_volcano.png", plot = soil_volcano, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "png")

ggsave("soil_volcano.pdf", plot = soil_volcano, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "pdf")


#For Rhizosphere sample

# Create -log10 pvalue safely
res_sig_rhizo_1$logP <- -log10(res_sig_rhizo_1$pvalue)

# Remove infinite values if any
res_sig_rhizo_1 <- res_sig_rhizo_1[is.finite(res_sig_rhizo_1$logP), ]

# Define differential expression category
res_sig_rhizo_1$diffexpressed <- "Not Significant"

res_sig_rhizo_1$diffexpressed[
  res_sig_rhizo_1$log2FoldChange > 0.6 & res_sig_rhizo_1$pvalue < 1e-10
] <- "GDD_1400"

res_sig_rhizo_1$diffexpressed[
  res_sig_rhizo_1$log2FoldChange < -0.6 & res_sig_rhizo_1$pvalue < 1e-10
] <- "GDD_600"

# Label ONLY highly significant taxa
res_sig_rhizo_1$delabel <- ifelse(
  res_sig_rhizo_1$diffexpressed != "Not Significant",
  res_sig_rhizo_1$Genus,
  NA
)

# Volcano plot
rhizo_volcano <- ggplot(res_sig_rhizo_1,
                        aes(x = log2FoldChange,
                            y = logP,
                            color = diffexpressed)) +
  
  geom_point(size = 15, alpha = 0.9) +
  
  geom_vline(xintercept = c(-0.6, 0.6),
             color = "#2C6BED",
             linewidth = 0.8) +
  
  geom_hline(yintercept = -log10(1e-7),
             color = "#2C6BED",
             linewidth = 0.8) +
  
  geom_text_repel(
    data = subset(res_sig_rhizo_1, !is.na(delabel)),
    aes(label = delabel),
    size = 4.5,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  
  scale_color_manual(
    values = c(
      "GDD_600" = "#2C6BED",      
      "Not Significant" = "grey70",
      "GDD_1400" = "#D62728"      
    ),
    name = "More Abundant In:"
  ) +
  
  labs(
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~italic(p)~value)
  ) +
  
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text  = element_text(size = 14),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank()
  )

rhizo_volcano

#save plot
ggsave("rhizo_volcano.png", plot = rhizo_volcano, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "png")

ggsave("rhizo_volcano.pdf", plot = rhizo_volcano, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "pdf")



#For root sample

# Create safe -log10 pvalue
res_sig_root_1$logP <- -log10(res_sig_root_1$pvalue)

# Remove infinite values (if any)
res_sig_root_1 <- res_sig_root_1[is.finite(res_sig_root_1$logP), ]

# Define stronger statistical separation (keeping 0.6)
res_sig_root_1$diffexpressed <- "Not Significant"

res_sig_root_1$diffexpressed[
  res_sig_root_1$log2FoldChange > 0.6 & res_sig_root_1$pvalue < 1e-11
] <- "GDD_1400"

res_sig_root_1$diffexpressed[
  res_sig_root_1$log2FoldChange < -0.6 & res_sig_root_1$pvalue < 1e-11
] <- "GDD_600"

# Label only highly significant taxa
res_sig_root_1$delabel <- ifelse(
  abs(res_sig_root_1$log2FoldChange) > 0.6 &
    res_sig_root_1$pvalue < 1e-11,
  res_sig_root_1$Genus,
  NA
)

# Volcano plot
root_volcano <- ggplot(res_sig_root_1,
                       aes(x = log2FoldChange,
                           y = logP,
                           color = diffexpressed)) +
  
  geom_point(size = 15, alpha = 0.9) +
  
  geom_vline(xintercept = c(-0.6, 0.6),
             color = "#2C6BED",
             linewidth = 0.8) +
  
  geom_hline(yintercept = -log10(1e-9),
             color = "#2C6BED",
             linewidth = 0.8) +
  
  geom_text_repel(
    data = subset(res_sig_root_1, !is.na(delabel)),
    aes(label = delabel),
    size = 4.5,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  
  scale_color_manual(
    values = c(
      "GDD_600" = "#2C6BED",
      "Not Significant" = "grey70",
      "GDD_1400" = "#D62728"
    ),
    name = "More Abundant In:"
  ) +
  
  labs(
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~italic(p)~value)
  ) +
  
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text  = element_text(size = 14),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank()
  )

root_volcano

#save plot
ggsave("root_volcano.png", plot = root_volcano, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "png")

ggsave("root_volcano.pdf", plot = root_volcano, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "pdf")


#For Stalk sample

# Create safe -log10 pvalue
res_sig_stalk_1$logP <- -log10(res_sig_stalk_1$pvalue)

# Remove infinite values if any
res_sig_stalk_1 <- res_sig_stalk_1[is.finite(res_sig_stalk_1$logP), ]

# Define stronger statistical separation (keeping 0.6)
res_sig_stalk_1$diffexpressed <- "Not Significant"

res_sig_stalk_1$diffexpressed[
  res_sig_stalk_1$log2FoldChange > 0.6 & res_sig_stalk_1$pvalue < 1e-12
] <- "GDD_1400"

res_sig_stalk_1$diffexpressed[
  res_sig_stalk_1$log2FoldChange < -0.6 & res_sig_stalk_1$pvalue < 1e-12
] <- "GDD_600"

# Label only highly significant taxa
res_sig_stalk_1$delabel <- ifelse(
  abs(res_sig_stalk_1$log2FoldChange) > 0.6 &
    res_sig_stalk_1$pvalue < 1e-12,
  res_sig_stalk_1$Genus,
  NA
)

# Volcano plot
stalk_volcano <- ggplot(res_sig_stalk_1,
                        aes(x = log2FoldChange,
                            y = logP,
                            color = diffexpressed)) +
  
  geom_point(size = 15, alpha = 0.9) +
  
  geom_vline(xintercept = c(-0.6, 0.6),
             color = "#2C6BED",
             linewidth = 0.8) +
  
  geom_hline(yintercept = -log10(1e-9),
             color = "#2C6BED",
             linewidth = 0.8) +
  
  geom_text_repel(
    data = subset(res_sig_stalk_1, !is.na(delabel)),
    aes(label = delabel),
    size = 4.5,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  
  scale_color_manual(
    values = c(
      "GDD_600" = "#2C6BED",
      "Not Significant" = "grey70",
      "GDD_1400" = "#D62728"
    ),
    name = "More Abundant In:"
  ) +
  
  labs(
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~italic(p)~value)
  ) +
  
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text  = element_text(size = 14),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank()
  )

stalk_volcano


#save plot
ggsave("stalk_volcano.png", plot = stalk_volcano, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width =15, height = 10, units = c("in"), device = "png")

ggsave("stalk_volcano.pdf", plot = stalk_volcano, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "pdf")