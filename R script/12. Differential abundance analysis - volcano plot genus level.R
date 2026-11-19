library(phyloseq)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)

set.seed(1234)

alpha <- 0.001
lfc_cutoff <- 0.6


#For Soil Sample

#-------------------------------------------------
# Read phyloseq object
#-------------------------------------------------
ps_soil <- readRDS(
  "C:/Aduragbemi/Microbiome/RDS/soil_bac_sub_no_bad_Filtered_final_non_normalized.rds"
)

#convert GDD to a factor
sample_data(ps_soil)$Growing_Degree_Days <- factor(
  sample_data(ps_soil)$Growing_Degree_Days, levels = c("1400", "600"))

#-------------------------------------------------
# Convert to DESeq2
#-------------------------------------------------
dds_soil <- phyloseq_to_deseq2(ps_soil, ~ Growing_Degree_Days)

# Recommended for sparse microbiome datasets
dds_soil <- estimateSizeFactors(dds_soil, type="poscounts")

#-------------------------------------------------
# Filter low abundance ASVs
#-------------------------------------------------
keep_soil <- rowSums(counts(dds_soil, normalized=TRUE) >= 2) >= 3
dds_soil <- dds_soil[keep_soil,]

#-------------------------------------------------
# Run DESeq2
#-------------------------------------------------
dds_soil <- DESeq(dds_soil)

# Positive log2FC = GDD600
res_soil <- results(
  dds_soil,
  contrast=c("Growing_Degree_Days","600","1400"),
  alpha=alpha
)

res_soil <- res_soil[order(res_soil$padj, na.last=NA),]

#-------------------------------------------------
# Merge taxonomy with ALL results
#-------------------------------------------------
tax_soil <- as.data.frame(
  tax_table(ps_soil)[rownames(res_soil),]
)

res_all_soil <- cbind(as.data.frame(res_soil),tax_soil)
res_all_soil <- res_all_soil[,!duplicated(colnames(res_all_soil))]

write.csv(
  res_all_soil,
  "C:/Aduragbemi/Manuscript/Review 8/Differential abundance/res_soil_all.csv",
  row.names=TRUE)

saveRDS(
  res_all_soil,
  "C:/Aduragbemi/Microbiome/RDS/res_soil_all.rds")

#-------------------------------------------------
# Significant ASVs only
#-------------------------------------------------
res_sig_soil_1 <- subset(
  res_all_soil,
  !is.na(padj) &
    padj < alpha
)

write.csv(
  res_sig_soil_1,
  "C:/Aduragbemi/Manuscript/Review 8/Differential abundance/res_sig_soil.csv",
  row.names=TRUE)

saveRDS(
  res_sig_soil_1,
  "C:/Aduragbemi/Microbiome/RDS/res_sig_soil_1.rds")


#-------------------------------------------------
# Volcano plot uses ALL ASVs
#-------------------------------------------------
volcano_soil <- subset(
  res_all_soil,
  !is.na(padj)
)

# Create -log10 adjusted p-value
volcano_soil$logP <- -log10(volcano_soil$padj)

# Define differential abundance status
volcano_soil$Status <- "Not Significant"

# Positive log2FC = enriched in GDD 600
volcano_soil$Status[
  volcano_soil$padj < alpha &
    volcano_soil$log2FoldChange > lfc_cutoff
] <- "GDD_600"

# Negative log2FC = enriched in GDD 1400
volcano_soil$Status[
  volcano_soil$padj < alpha &
    volcano_soil$log2FoldChange < -lfc_cutoff
] <- "GDD_1400"

#-------------------------------------------------
# Label ONLY the 5 most significant taxa from each side
#-------------------------------------------------

volcano_soil$Label <- NA

top600 <- volcano_soil %>%
  filter(Status == "GDD_600") %>%
  arrange(padj) %>%
  slice_head(n = 5)

top1400 <- volcano_soil %>%
  filter(Status == "GDD_1400") %>%
  arrange(padj) %>%
  slice_head(n = 5)

top_labels <- bind_rows(top600, top1400)

volcano_soil$Label[
  rownames(volcano_soil) %in% rownames(top_labels)
] <- volcano_soil$Genus[
  rownames(volcano_soil) %in% rownames(top_labels)
]

#-------------------------------------------------
# Volcano plot
#-------------------------------------------------

soil_volcano <-
  ggplot(
    volcano_soil,
    aes(
      x = log2FoldChange,
      y = logP,
      colour = Status
    )
  ) +
  
  geom_point(
    size = 15,
    alpha = 0.85
  ) +
  
  geom_vline(
    xintercept = c(-lfc_cutoff, lfc_cutoff),
    linetype = "dashed",
    linewidth = 0.8,
    colour = "#2C6BED"
  ) +
  
  geom_hline(
    yintercept = -log10(alpha),
    linetype = "dashed",
    linewidth = 0.8,
    colour = "#2C6BED"
  ) +
  
  geom_text_repel(
    data = subset(volcano_soil, !is.na(Label)),
    aes(label = Label),
    size = 17,
    box.padding = 0.5,
    point.padding = 0.4,
    segment.size = 0.3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  
  scale_colour_manual(
    values = c(
      "GDD_600" = "#1F78B4",
      "GDD_1400" = "#E31A1C",
      "Not Significant" = "grey75"
    ),
    name = "More abundant in"
  ) +
  
  labs(
    title = "Soil",
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~Adjusted~italic(P))
  ) +
  
  theme_bw(base_size = 16) +
  
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 25),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

soil_volcano

#-------------------------------------------------
# Save figure
#-------------------------------------------------

ggsave(
  "C:/Aduragbemi/Manuscript/Review 8/IMAGES8/soil_volcano.png",
  soil_volcano,
  width = 15,
  height = 10,
  dpi = 700
)

ggsave(
  "C:/Aduragbemi/Manuscript/Review 8/PDF8/soil_volcano.pdf",
  soil_volcano,
  width = 15,
  height = 10,
  dpi = 700
)


#For Rhizosphere Sample

#-------------------------------------------------
# Read phyloseq object
#-------------------------------------------------
ps_rhizo <- readRDS(
  "C:/Aduragbemi/Microbiome/RDS/rhizo_bac_sub_no_bad_Filtered_final_non_normalized.rds"
)

#convert GDD to a factor
sample_data(ps_rhizo)$Growing_Degree_Days <- factor(
  sample_data(ps_rhizo)$Growing_Degree_Days, levels = c("1400", "600"))

#-------------------------------------------------
# Convert to DESeq2
#-------------------------------------------------
dds_rhizo <- phyloseq_to_deseq2(ps_rhizo, ~ Growing_Degree_Days)

# Recommended for sparse microbiome datasets
dds_rhizo <- estimateSizeFactors(dds_rhizo, type="poscounts")

#-------------------------------------------------
# Filter low abundance ASVs
#-------------------------------------------------
keep_rhizo <- rowSums(counts(dds_rhizo, normalized=TRUE) >= 2) >= 3
dds_rhizo <- dds_rhizo[keep_rhizo,]

#-------------------------------------------------
# Run DESeq2
#-------------------------------------------------
dds_rhizo <- DESeq(dds_rhizo)

# Positive log2FC = GDD600
res_rhizo <- results(
  dds_rhizo,
  contrast=c("Growing_Degree_Days","600","1400"),
  alpha=alpha
)

res_rhizo <- res_rhizo[order(res_rhizo$padj, na.last=NA),]

#-------------------------------------------------
# Merge taxonomy with ALL results
#-------------------------------------------------
tax_rhizo <- as.data.frame(
  tax_table(ps_rhizo)[rownames(res_rhizo),]
)

res_all_rhizo <- cbind(as.data.frame(res_rhizo),tax_rhizo)
res_all_rhizo <- res_all_rhizo[,!duplicated(colnames(res_all_rhizo))]

write.csv(
  res_all_rhizo,
  "C:/Aduragbemi/Manuscript/Review 8/Differential abundance/res_rhizo_all.csv",
  row.names=TRUE)

saveRDS(
  res_all_rhizo,
  "C:/Aduragbemi/Microbiome/RDS/res_rhizo_all.rds")

#-------------------------------------------------
# Significant ASVs only
#-------------------------------------------------
res_sig_rhizo_1 <- subset(
  res_all_rhizo,
  !is.na(padj) &
    padj < alpha
)

write.csv(
  res_sig_rhizo_1,
  "C:/Aduragbemi/Manuscript/Review 8/Differential abundance/res_sig_rhizo.csv",
  row.names=TRUE)

saveRDS(
  res_sig_rhizo_1,
  "C:/Aduragbemi/Microbiome/RDS/res_sig_rhizo_1.rds")


#-------------------------------------------------
# Volcano plot uses ALL ASVs
#-------------------------------------------------
volcano_rhizo <- subset(
  res_all_rhizo,
  !is.na(padj)
)

# Create -log10 adjusted p-value
volcano_rhizo$logP <- -log10(volcano_rhizo$padj)

# Define differential abundance status
volcano_rhizo$Status <- "Not Significant"

# Positive log2FC = enriched in GDD 600
volcano_rhizo$Status[
  volcano_rhizo$padj < alpha &
    volcano_rhizo$log2FoldChange > lfc_cutoff
] <- "GDD_600"

# Negative log2FC = enriched in GDD 1400
volcano_rhizo$Status[
  volcano_rhizo$padj < alpha &
    volcano_rhizo$log2FoldChange < -lfc_cutoff
] <- "GDD_1400"

#-------------------------------------------------
# Label ONLY the 5 most significant taxa from each side
#-------------------------------------------------

volcano_rhizo$Label <- NA

top600 <- volcano_rhizo %>%
  filter(Status == "GDD_600") %>%
  arrange(padj) %>%
  slice_head(n = 5)

top1400 <- volcano_rhizo %>%
  filter(Status == "GDD_1400") %>%
  arrange(padj) %>%
  slice_head(n = 5)

top_labels <- bind_rows(top600, top1400)

volcano_rhizo$Label[
  rownames(volcano_rhizo) %in% rownames(top_labels)
] <- volcano_rhizo$Genus[
  rownames(volcano_rhizo) %in% rownames(top_labels)
]

#-------------------------------------------------
# Volcano plot
#-------------------------------------------------

rhizo_volcano <-
  ggplot(
    volcano_rhizo,
    aes(
      x = log2FoldChange,
      y = logP,
      colour = Status
    )
  ) +
  
  geom_point(
    size = 15,
    alpha = 0.85
  ) +
  
  geom_vline(
    xintercept = c(-lfc_cutoff, lfc_cutoff),
    linetype = "dashed",
    linewidth = 0.8,
    colour = "#2C6BED"
  ) +
  
  geom_hline(
    yintercept = -log10(alpha),
    linetype = "dashed",
    linewidth = 0.8,
    colour = "#2C6BED"
  ) +
  
  geom_text_repel(
    data = subset(volcano_rhizo, !is.na(Label)),
    aes(label = Label),
    size = 17,
    box.padding = 0.5,
    point.padding = 0.4,
    segment.size = 0.3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  
  scale_colour_manual(
    values = c(
      "GDD_600" = "#1F78B4",
      "GDD_1400" = "#E31A1C",
      "Not Significant" = "grey75"
    ),
    name = "More abundant in"
  ) +
  
  labs(
    title = "Rhizosphere",
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~Adjusted~italic(P))
  ) +
  
  theme_bw(base_size = 16) +
  
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 25),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

rhizo_volcano



#-------------------------------------------------
# Save figure
#-------------------------------------------------

ggsave(
  "C:/Aduragbemi/Manuscript/Review 8/IMAGES8/rhizo_volcano.png",
  rhizo_volcano,
  width = 15,
  height = 10,
  dpi = 700
)

ggsave(
  "C:/Aduragbemi/Manuscript/Review 8/PDF8/rhizo_volcano.pdf",
  rhizo_volcano,
  width = 15,
  height = 10,
  dpi = 700
)



#For Root Sample

#-------------------------------------------------
# Read phyloseq object
#-------------------------------------------------
ps_root <- readRDS(
  "C:/Aduragbemi/Microbiome/RDS/root_bac_sub_no_bad_Filtered_final_non_normalized.rds"
)

#convert GDD to a factor
sample_data(ps_root)$Growing_Degree_Days <- factor(
  sample_data(ps_root)$Growing_Degree_Days, levels = c("1400", "600"))

#-------------------------------------------------
# Convert to DESeq2
#-------------------------------------------------
dds_root <- phyloseq_to_deseq2(ps_root, ~ Growing_Degree_Days)

# Recommended for sparse microbiome datasets
dds_root <- estimateSizeFactors(dds_root, type="poscounts")

#-------------------------------------------------
# Filter low abundance ASVs
#-------------------------------------------------
keep_root <- rowSums(counts(dds_root, normalized=TRUE) >= 2) >= 3
dds_root <- dds_root[keep_root,]

#-------------------------------------------------
# Run DESeq2
#-------------------------------------------------
dds_root <- DESeq(dds_root)

# Positive log2FC = GDD600
res_root <- results(
  dds_root,
  contrast=c("Growing_Degree_Days","600","1400"),
  alpha=alpha
)

res_root <- res_root[order(res_root$padj, na.last=NA),]

#-------------------------------------------------
# Merge taxonomy with ALL results
#-------------------------------------------------
tax_root <- as.data.frame(
  tax_table(ps_root)[rownames(res_root),]
)

res_all_root <- cbind(as.data.frame(res_root),tax_root)
res_all_root <- res_all_root[,!duplicated(colnames(res_all_root))]

write.csv(
  res_all_root,
  "C:/Aduragbemi/Manuscript/Review 8/Differential abundance/res_root_all.csv",
  row.names=TRUE)

saveRDS(
  res_all_root,
  "C:/Aduragbemi/Microbiome/RDS/res_root_all.rds")

#-------------------------------------------------
# Significant ASVs only
#-------------------------------------------------
res_sig_root_1 <- subset(
  res_all_root,
  !is.na(padj) &
    padj < alpha
)

write.csv(
  res_sig_root_1,
  "C:/Aduragbemi/Manuscript/Review 8/Differential abundance/res_sig_root.csv",
  row.names=TRUE)

saveRDS(
  res_sig_root_1,
  "C:/Aduragbemi/Microbiome/RDS/res_sig_root_1.rds")


#-------------------------------------------------
# Volcano plot uses ALL ASVs
#-------------------------------------------------
volcano_root <- subset(
  res_all_root,
  !is.na(padj)
)

# Create -log10 adjusted p-value
volcano_root$logP <- -log10(volcano_root$padj)

# Define differential abundance status
volcano_root$Status <- "Not Significant"

# Positive log2FC = enriched in GDD 600
volcano_root$Status[
  volcano_root$padj < alpha &
    volcano_root$log2FoldChange > lfc_cutoff
] <- "GDD_600"

# Negative log2FC = enriched in GDD 1400
volcano_root$Status[
  volcano_root$padj < alpha &
    volcano_root$log2FoldChange < -lfc_cutoff
] <- "GDD_1400"

#-------------------------------------------------
# Label ONLY the 5 most significant taxa from each side
#-------------------------------------------------

volcano_root$Label <- NA

top600 <- volcano_root %>%
  filter(Status == "GDD_600") %>%
  arrange(padj) %>%
  slice_head(n = 5)

top1400 <- volcano_root %>%
  filter(Status == "GDD_1400") %>%
  arrange(padj) %>%
  slice_head(n = 5)

top_labels <- bind_rows(top600, top1400)

volcano_root$Label[
  rownames(volcano_root) %in% rownames(top_labels)
] <- volcano_root$Genus[
  rownames(volcano_root) %in% rownames(top_labels)
]

#-------------------------------------------------
# Volcano plot
#-------------------------------------------------

root_volcano <-
  ggplot(
    volcano_root,
    aes(
      x = log2FoldChange,
      y = logP,
      colour = Status
    )
  ) +
  
  geom_point(
    size = 15,
    alpha = 0.85
  ) +
  
  geom_vline(
    xintercept = c(-lfc_cutoff, lfc_cutoff),
    linetype = "dashed",
    linewidth = 0.8,
    colour = "#2C6BED"
  ) +
  
  geom_hline(
    yintercept = -log10(alpha),
    linetype = "dashed",
    linewidth = 0.8,
    colour = "#2C6BED"
  ) +
  
  geom_text_repel(
    data = subset(volcano_root, !is.na(Label)),
    aes(label = Label),
    size = 17,
    box.padding = 0.5,
    point.padding = 0.4,
    segment.size = 0.3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  
  scale_colour_manual(
    values = c(
      "GDD_600" = "#1F78B4",
      "GDD_1400" = "#E31A1C",
      "Not Significant" = "grey75"
    ),
    name = "More abundant in"
  ) +
  
  labs(
    title = "root",
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~Adjusted~italic(P))
  ) +
  
  theme_bw(base_size = 16) +
  
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 25),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

root_volcano



#-------------------------------------------------
# Save figure
#-------------------------------------------------

ggsave(
  "C:/Aduragbemi/Manuscript/Review 8/IMAGES8/root_volcano.png",
  root_volcano,
  width = 15,
  height = 10,
  dpi = 700
)

ggsave(
  "C:/Aduragbemi/Manuscript/Review 8/PDF8/root_volcano.pdf",
  root_volcano,
  width = 15,
  height = 10,
  dpi = 700
)





#For Stalk Sample

#-------------------------------------------------
# Read phyloseq object
#-------------------------------------------------
ps_stalk <- readRDS(
  "C:/Aduragbemi/Microbiome/RDS/stalk_bac_sub_no_bad_Filtered_final_non_normalized.rds"
)

#convert GDD to a factor
sample_data(ps_stalk)$Growing_Degree_Days <- factor(
  sample_data(ps_stalk)$Growing_Degree_Days, levels = c("1400", "600"))

#-------------------------------------------------
# Convert to DESeq2
#-------------------------------------------------
dds_stalk <- phyloseq_to_deseq2(ps_stalk, ~ Growing_Degree_Days)

# Recommended for sparse microbiome datasets
dds_stalk <- estimateSizeFactors(dds_stalk, type="poscounts")

#-------------------------------------------------
# Filter low abundance ASVs
#-------------------------------------------------
keep_stalk <- rowSums(counts(dds_stalk, normalized=TRUE) >= 2) >= 3
dds_stalk <- dds_stalk[keep_stalk,]

#-------------------------------------------------
# Run DESeq2
#-------------------------------------------------
dds_stalk <- DESeq(dds_stalk)

# Positive log2FC = GDD600
res_stalk <- results(
  dds_stalk,
  contrast=c("Growing_Degree_Days","600","1400"),
  alpha=alpha
)

res_stalk <- res_stalk[order(res_stalk$padj, na.last=NA),]

#-------------------------------------------------
# Merge taxonomy with ALL results
#-------------------------------------------------
tax_stalk <- as.data.frame(
  tax_table(ps_stalk)[rownames(res_stalk),]
)

res_all_stalk <- cbind(as.data.frame(res_stalk),tax_stalk)
res_all_stalk <- res_all_stalk[,!duplicated(colnames(res_all_stalk))]

write.csv(
  res_all_stalk,
  "C:/Aduragbemi/Manuscript/Review 8/Differential abundance/res_stalk_all.csv",
  row.names=TRUE)

saveRDS(
  res_all_stalk,
  "C:/Aduragbemi/Microbiome/RDS/res_stalk_all.rds")

#-------------------------------------------------
# Significant ASVs only
#-------------------------------------------------
res_sig_stalk_1 <- subset(
  res_all_stalk,
  !is.na(padj) &
    padj < alpha
)

write.csv(
  res_sig_stalk_1,
  "C:/Aduragbemi/Manuscript/Review 8/Differential abundance/res_sig_stalk.csv",
  row.names=TRUE)

saveRDS(
  res_sig_stalk_1,
  "C:/Aduragbemi/Microbiome/RDS/res_sig_stalk_1.rds")


#-------------------------------------------------
# Volcano plot uses ALL ASVs
#-------------------------------------------------
volcano_stalk <- subset(
  res_all_stalk,
  !is.na(padj)
)

# Create -log10 adjusted p-value
volcano_stalk$logP <- -log10(volcano_stalk$padj)

# Define differential abundance status
volcano_stalk$Status <- "Not Significant"

# Positive log2FC = enriched in GDD 600
volcano_stalk$Status[
  volcano_stalk$padj < alpha &
    volcano_stalk$log2FoldChange > lfc_cutoff
] <- "GDD_600"

# Negative log2FC = enriched in GDD 1400
volcano_stalk$Status[
  volcano_stalk$padj < alpha &
    volcano_stalk$log2FoldChange < -lfc_cutoff
] <- "GDD_1400"

#-------------------------------------------------
# Label ONLY the 5 most significant taxa from each side
#-------------------------------------------------

volcano_stalk$Label <- NA

top600 <- volcano_stalk %>%
  filter(Status == "GDD_600") %>%
  arrange(padj) %>%
  slice_head(n = 5)

top1400 <- volcano_stalk %>%
  filter(Status == "GDD_1400") %>%
  arrange(padj) %>%
  slice_head(n = 5)

top_labels <- bind_rows(top600, top1400)

volcano_stalk$Label[
  rownames(volcano_stalk) %in% rownames(top_labels)
] <- volcano_stalk$Genus[
  rownames(volcano_stalk) %in% rownames(top_labels)
]

#-------------------------------------------------
# Volcano plot
#-------------------------------------------------

stalk_volcano <-
  ggplot(
    volcano_stalk,
    aes(
      x = log2FoldChange,
      y = logP,
      colour = Status
    )
  ) +
  
  geom_point(
    size = 15,
    alpha = 0.85
  ) +
  
  geom_vline(
    xintercept = c(-lfc_cutoff, lfc_cutoff),
    linetype = "dashed",
    linewidth = 0.8,
    colour = "#2C6BED"
  ) +
  
  geom_hline(
    yintercept = -log10(alpha),
    linetype = "dashed",
    linewidth = 0.8,
    colour = "#2C6BED"
  ) +
  
  geom_text_repel(
    data = subset(volcano_stalk, !is.na(Label)),
    aes(label = Label),
    size = 17,
    box.padding = 0.5,
    point.padding = 0.4,
    segment.size = 0.3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  
  scale_colour_manual(
    values = c(
      "GDD_600" = "#1F78B4",
      "GDD_1400" = "#E31A1C",
      "Not Significant" = "grey75"
    ),
    name = "More abundant in"
  ) +
  
  labs(
    title = "stalk",
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~Adjusted~italic(P))
  ) +
  
  theme_bw(base_size = 16) +
  
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 25),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

stalk_volcano



#-------------------------------------------------
# Save figure
#-------------------------------------------------

ggsave(
  "C:/Aduragbemi/Manuscript/Review 8/IMAGES8/stalk_volcano.png",
  stalk_volcano,
  width = 15,
  height = 10,
  dpi = 700
)

ggsave(
  "C:/Aduragbemi/Manuscript/Review 8/PDF8/stalk_volcano.pdf",
  stalk_volcano,
  width = 15,
  height = 10,
  dpi = 700
)
