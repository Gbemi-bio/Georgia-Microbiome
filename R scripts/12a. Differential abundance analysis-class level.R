# Load required libraries
library(phyloseq)
library(DESeq2)
library(ggplot2)

set.seed(1234)


# ================================
# read in the non-normalized data
# ================================
root_bac_sub_no_bad_Filtered_final_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/root_bac_sub_no_bad_Filtered_final_non_normalized.rds")
stalk_bac_sub_no_bad_Filtered_final_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/stalk_bac_sub_no_bad_Filtered_final_non_normalized.rds")
soil_bac_sub_no_bad_Filtered_final_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/soil_bac_sub_no_bad_Filtered_final_non_normalized.rds")
rhizo_bac_sub_no_bad_Filtered_final_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/rhizo_bac_sub_no_bad_Filtered_final_non_normalized.rds")

root_bac_sub_no_bad_Filtered_final_non_normalized
stalk_bac_sub_no_bad_Filtered_final_non_normalized 
soil_bac_sub_no_bad_Filtered_final_non_normalized 
rhizo_bac_sub_no_bad_Filtered_final_non_normalized


# =========================================================
# Differential abundance analysis of sample types by class
# =========================================================

#For Soil sample

#apply a transformation to every sample in my ASV table
diff_soil <- transform_sample_counts(soil_bac_sub_no_bad_Filtered_final_non_normalized, function(OTU) OTU +1)

#Conversion to DESeq2 object
diff_soil = phyloseq_to_deseq2(diff_soil, ~ Growing_Degree_Days)

# calculate geometric means prior to estimate size factors
gm_mean_soil = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 1.0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diff_soil), 1, gm_mean_soil)

#filter out low count ASV before running DESeq, filter out ASVs where there are less than 3 samples with normalized counts greater than or equal to 3
diff_soil_est <- estimateSizeFactors(diff_soil, geoMeans = geoMeans)
diff_soil_est_filter <- rowSums( counts(diff_soil_est, normalized=TRUE) >= 2 ) >= 3
ds_soil <- diff_soil[diff_soil_est_filter,]
ds_soil

#Negative Binomial log-linear model with DESeq2
ds_soil = DESeq(ds_soil )
ds_soil
alpha = 0.001
res_soil = results(ds_soil, contrast=c("Growing_Degree_Days", "600", "1400"), alpha=alpha)
res_soil
res_soil = res_soil[order(res_soil$padj, na.last=NA), ]
res_sig_soil = res_soil[(res_soil$padj < alpha), ]
res_sig_soil
res_sig_soil = cbind(as(res_sig_soil, "data.frame"), as(phyloseq::tax_table(soil_bac_sub_no_bad_Filtered_final_non_normalized)[rownames(res_sig_soil), ], "matrix"))
res_sig_soil
head(res_sig_soil)

#remove duplicated columns
res_sig_soil_1 <- res_sig_soil[ , !duplicated(colnames(res_sig_soil))] 
res_sig_soil_1

#write out the result in res_sig_soil_1 to a cvs file
write.csv(as.data.frame(res_sig_soil_1), file="C:/Aduragbemi/Manuscript/Review 5/Differential abundance/res_sig_soil.csv")

#Save as RDS
saveRDS(res_sig_soil_1, file = "C:/Aduragbemi/Microbiome/RDS/res_sig_soil_1.rds")

#Plot the graph
theme_set(theme_grey())
sigtabgen_soil = subset(res_sig_soil_1, !is.na(Genus))

final_soil_plot <- ggplot(res_sig_soil_1) + 
  geom_bar(aes(x = log2FoldChange, y = Class, fill = ifelse(log2FoldChange > 0, "tomato2", "darkcyan")), 
           width = 1, 
           stat = "identity", 
           color = "black", 
           size = 1.5) +  # Adjust the size parameter for thicker outlines
  theme(axis.title = element_text(family = "serif", size = 40, colour = "black", face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 40, colour = "black"),
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.3, size = 40, colour = "black"),
        legend.title = element_text(colour = "black", size = 35, face = "bold"),
        legend.text = element_text(size = 35, colour = "black")) +
  geom_vline(xintercept = 0.0, color = "black", size = 1) +
  scale_fill_identity(guide = "none") +  
  guides(color = guide_legend(override.aes = list(size = 3)))

final_soil_plot

#save plot
ggsave("final_soil_plot.png", plot = final_soil_plot, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 30, height = 25, units = c("in"), device = "png")

ggsave("final_soil_plot.pdf", plot = final_soil_plot, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 30, height = 25, units = c("in"), device = "pdf")


#for Rhizosphere sample

#apply a transformation to every sample in my ASV table
diff_rhizo <- transform_sample_counts(rhizo_bac_sub_no_bad_Filtered_final_non_normalized, function(OTU) OTU +1)

#Conversion to DESeq2 object
diff_rhizo = phyloseq_to_deseq2(diff_rhizo, ~ Growing_Degree_Days)

# calculate geometric means prior to estimate size factors
gm_mean_rhizo = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 1.0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diff_rhizo), 1, gm_mean_rhizo)

#filter out low count ASV before running DESeq, filter out ASVs where there are less than 3 samples with normalized counts greater than or equal to 3
diff_rhizo_est <- estimateSizeFactors(diff_rhizo, geoMeans = geoMeans)
diff_rhizo_est_filter <- rowSums( counts(diff_rhizo_est, normalized=TRUE) >= 2 ) >= 3
ds_rhizo <- diff_rhizo[diff_rhizo_est_filter,]
ds_rhizo

#Negative Binomial log-linear model with DESeq2
ds_rhizo = DESeq(ds_rhizo )
ds_rhizo
alpha = 0.001
res_rhizo = results(ds_rhizo, contrast=c("Growing_Degree_Days", "600", "1400"), alpha=alpha)
res_rhizo
res_rhizo = res_rhizo[order(res_rhizo$padj, na.last=NA), ]
res_sig_rhizo = res_rhizo[(res_rhizo$padj < alpha), ]
res_sig_rhizo
res_sig_rhizo = cbind(as(res_sig_rhizo, "data.frame"), as(phyloseq::tax_table(rhizo_bac_sub_no_bad_Filtered_final_non_normalized)[rownames(res_sig_rhizo), ], "matrix"))
res_sig_rhizo

#remove duplicated columns
res_sig_rhizo_1 <- res_sig_rhizo[ , !duplicated(colnames(res_sig_rhizo))] 
res_sig_rhizo_1

#write it out the result in res_sig_rhizo_1 out in csv
write.csv(as.data.frame(res_sig_rhizo_1), file="C:/Aduragbemi/Manuscript/Review 5/Differential abundance/res_sig_rhizo.csv")

#Save as RDS
saveRDS(res_sig_rhizo_1, file = "C:/Aduragbemi/Microbiome/RDS/res_sig_rhizo_1.rds")

#plotting graph
theme_set(theme_grey())
sigtabgen_rhizo = subset(res_sig_rhizo_1, !is.na(Class))
final_rhizo_plot <- ggplot(res_sig_rhizo_1) + 
  geom_bar(aes(x = log2FoldChange, y = Class, fill = ifelse(log2FoldChange > 0, "tomato2", "darkcyan")), 
           width = 1, 
           stat = "identity", 
           color = "black", 
           size = 1.5) +  # Adjust the size parameter for thicker outlines
  theme(axis.title = element_text(family = "serif", size = 40, colour = "black", face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 40, colour = "black"),
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.3, size = 40, colour = "black"),
        legend.title = element_text(colour = "black", size = 35, face = "bold"),
        legend.text = element_text(size = 35, colour = "black")) +
  geom_vline(xintercept = 0.0, color = "black", size = 1) +
  scale_fill_identity(guide = "none") + 
  guides(color = guide_legend(override.aes = list(size = 3)))

final_rhizo_plot

#save plot
ggsave("final_rhizo_plot.png", plot = final_rhizo_plot, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 30, height = 25, units = c("in"), device = "png")

ggsave("final_rhizo_plot.pdf", plot = final_rhizo_plot, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 30, height = 25, units = c("in"), device = "pdf")


#For Root sample

#apply a transformation to every sample in my ASV table
diff_root <- transform_sample_counts(root_bac_sub_no_bad_Filtered_final_non_normalized, function(OTU) OTU +1)

#Conversion to DESeq2 object
diff_root = phyloseq_to_deseq2(diff_root, ~ Growing_Degree_Days)

# calculate geometric means prior to estimate size factors
gm_mean_root = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 1.0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diff_root), 1, gm_mean_root)

#filter out low count ASV before running DESeq, filter out ASVs where there are less than 3 samples with normalized counts greater than or equal to 3
diff_root_est <- estimateSizeFactors(diff_root, geoMeans = geoMeans)
diff_root_est_filter <- rowSums( counts(diff_root_est, normalized=TRUE) >= 2 ) >= 3
ds_root <- diff_root[diff_root_est_filter,]
ds_root

#Negative Binomial log-linear model with DESeq2
ds_root = DESeq(ds_root )
ds_root
alpha = 0.001
res_root = results(ds_root, contrast=c("Growing_Degree_Days", "600", "1400"), alpha=alpha)
res_root
res_root = res_root[order(res_root$padj, na.last=NA), ]
res_sig_root = res_root[(res_root$padj < alpha), ]
res_sig_root
res_sig_root = cbind(as(res_sig_root, "data.frame"), as(phyloseq::tax_table(root_bac_sub_no_bad_Filtered_final_non_normalized)[rownames(res_sig_root), ], "matrix"))
res_sig_root

#remove duplicated columns
res_sig_root_1 <- res_sig_root[ , !duplicated(colnames(res_sig_root))] 
res_sig_root_1

#write it out the result in res_sig_root_1 out in csv
write.csv(as.data.frame(res_sig_root_1), file="C:/Aduragbemi/Manuscript/Review 5/Differential abundance/res_sig_root.csv")

#Save as RDS
saveRDS(res_sig_root_1, file = "C:/Aduragbemi/Microbiome/RDS/res_sig_root_1.rds")

#plotting graph
theme_set(theme_grey())
sigtabgen_root = subset(res_sig_root_1, !is.na(Class))
final_root_plot <- ggplot(res_sig_root_1) + 
  geom_bar(aes(x = log2FoldChange, y = Class, fill = ifelse(log2FoldChange > 0, "tomato2", "darkcyan")), 
           width = 1, 
           stat = "identity", 
           color = "black", 
           size = 1.5) +  # Adjust the size parameter for thicker outlines
  theme(axis.title = element_text(family = "serif", size = 40, colour = "black", face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 40, colour = "black"),
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.3, size = 40, colour = "black"),
        legend.title = element_text(colour = "black", size = 35, face = "bold"),
        legend.text = element_text(size = 35, colour = "black")) +
  geom_vline(xintercept = 0.0, color = "black", size = 1) +
  scale_fill_identity(guide = "none") + 
  guides(color = guide_legend(override.aes = list(size = 3)))

final_root_plot

#save plot
ggsave("final_root_plot.png", plot = final_root_plot, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 30, height = 25, units = c("in"), device = "png")

ggsave("final_root_plot.pdf", plot = final_root_plot, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 30, height = 25, units = c("in"), device = "pdf")



#For Stalk sample

#apply a transformation to every sample in my ASV table
diff_stalk <- transform_sample_counts(stalk_bac_sub_no_bad_Filtered_final_non_normalized, function(OTU) OTU +1)

#Conversion to DESeq2 object
diff_stalk = phyloseq_to_deseq2(diff_stalk, ~ Growing_Degree_Days)

# calculate geometric means prior to estimate size factors
gm_mean_stalk = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 1.0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diff_stalk), 1, gm_mean_stalk)

#filter out low count ASV before running DESeq, filter out ASVs where there are less than 3 samples with normalized counts greater than or equal to 3
diff_stalk_est <- estimateSizeFactors(diff_stalk, geoMeans = geoMeans)
diff_stalk_est_filter <- rowSums( counts(diff_stalk_est, normalized=TRUE) >= 2 ) >= 3
ds_stalk <- diff_stalk[diff_stalk_est_filter,]
ds_stalk

#Negative Binomial log-linear model with DESeq2
ds_stalk = DESeq(ds_stalk )
ds_stalk
alpha = 0.001
res_stalk = results(ds_stalk, contrast=c("Growing_Degree_Days", "600", "1400"), alpha=alpha)
res_stalk
res_stalk = res_stalk[order(res_stalk$padj, na.last=NA), ]
res_sig_stalk = res_stalk[(res_stalk$padj < alpha), ]
res_sig_stalk
res_sig_stalk = cbind(as(res_sig_stalk, "data.frame"), as(phyloseq::tax_table(stalk_bac_sub_no_bad_Filtered_final_non_normalized)[rownames(res_sig_stalk), ], "matrix"))
res_sig_stalk

#remove duplicated columns
res_sig_stalk_1 <- res_sig_stalk[ , !duplicated(colnames(res_sig_stalk))] 
res_sig_stalk_1

#write it out the result in res_sig_stalk_1 out in csv
write.csv(as.data.frame(res_sig_stalk_1), file="C:/Aduragbemi
          /Manuscript/Review 5/Differential abundance/res_sig_stalk.csv")

#Save as RDS
saveRDS(res_sig_stalk_1, file = "C:/Aduragbemi/Microbiome/RDS/res_sig_stalk_1.rds")


#plotting graph
theme_set(theme_grey())
sigtabgen_stalk = subset(res_sig_stalk_1, !is.na(Class))
final_stalk_plot <- ggplot(res_sig_stalk_1) + 
  geom_bar(aes(x = log2FoldChange, y = Class, fill = ifelse(log2FoldChange > 0, "tomato2", "darkcyan")), 
           width = 1, 
           stat = "identity", 
           color = "black", 
           size = 1.5) +  # Adjust the size parameter for thicker outlines
  theme(axis.title = element_text(family = "serif", size = 40, colour = "black", face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 40, colour = "black"),
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.3, size = 40, colour = "black"),
        legend.title = element_text(colour = "black", size = 35, face = "bold"),
        legend.text = element_text(size = 35, colour = "black")) +
  geom_vline(xintercept = 0.0, color = "black", size = 1) +
  scale_fill_identity(guide = "none") + 
  guides(color = guide_legend(override.aes = list(size = 3)))

final_stalk_plot

#save plot
ggsave("final_stalk_plot.png", plot = final_stalk_plot, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 30, height = 25, units = c("in"), device = "png")

ggsave("final_stalk_plot.pdf", plot = final_stalk_plot, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 30, height = 25, units = c("in"), device = "pdf")
