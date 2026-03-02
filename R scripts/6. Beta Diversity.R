#Load Libraries
library(phyloseq)
library(rbiom)
library(vegan)
library(ggplot2)
library(ggpubr)
library(ape)

# =============
# Define colors
# =============

colors_1 <- c("#008080", "#D2691E", "#6A5ACD", "#FF1493")
colors_2 <- c("#7FCDBB", "#006D6F")

# ==========================================================================
# Differences in microbial composition (beta-diversity) for all sample types
# ==========================================================================

#Read in the non normalized OTUs for all sample type together here
all_samples_non_normalized_final <- readRDS(file = "C:/Aduragbemi/Microbiome/RDs/all_samples_non_normalized_final.rds")
ps <- all_samples_non_normalized_final
ps_tree = phy_tree(ps)

# Check if binary (dichotomy) & multifurcating (polytomy) trees
sprintf("Is tree binary: %s", is.binary(ps_tree))

# If FASLE, randomly resolve polytomies and replace tree in "ps"
phy_tree(ps) = multi2di(ps_tree)
sprintf("Is tree binary: %s", is.binary(phy_tree(ps)))

#reorder the facet grid 
ps@sam_data[["Sample_Type"]] <- factor(ps@sam_data[["Sample_Type"]], 
                                       levels=c("Soil", "Rhizosphere", "Root", "Stalk"))

# Convert phyloseq object to rbiom object
biom_all <- as_rbiom(ps)
biom_all

#view th biom file
glimpse(biom_all, width = NULL)

#calculate the unifrac distance between samples
biom_dm <- bdiv_distmat(biom_all, 'unifrac')
biom_dm

#run a Principal Coordinates Analysis on my distance matrix.
pcoa_biom <- ape::pcoa(biom_dm)
pcoa_biom

#Calculate the Percent variance
PCoA1_biom <- round(pcoa_biom$values$Relative_eig[1] * 100, 2)
PCoA2_biom <- round(pcoa_biom$values$Relative_eig[2] * 100, 2)

#Ordinate samples based on beta diversity distances
wt.unifrac_all <- bdiv_ord_plot(
  biom_all,
  layers   = "pemt",
  bdiv     = "UniFrac",
  stat.by  = "Sample_Type",
  rank     = NULL
)

# Extract the ordination data
ord_data <- wt.unifrac_all$data

# Check columns
colnames(ord_data)

# Add points and title
wt.unifrac_all <- wt.unifrac_all +
  guides(shape = "none") +
  ggtitle("Weighted UniFrac Diversity PCoA") +
  geom_point(aes(color = Sample_Type), size = 6)

# Apply theme, legend, and color scale
wt.unifrac_all <- wt.unifrac_all +
  theme_gray() +
  theme(
    plot.title = element_text(family = "serif", face = "bold", size = 15),
    axis.title = element_text(family = "serif", size = 15, colour = "black", face = "bold"),
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 15, colour = "black"),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 15, colour = "black"),
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(size = 15, colour = "black"),
    legend.key.width = unit(1, "cm")) +
    xlab("PCoA 1 (30.7%)") +
    ylab("PCoA 2 (11.7%)") +
  scale_color_manual(values = colors_1, name = "Sample Types")

# Draw plot
wt.unifrac_all

#save plot as image
ggsave("wt.unifrac_all.png", plot = wt.unifrac_all , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save plot as pdf
ggsave("wt.unifrac_all.pdf", plot = wt.unifrac_all, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")



# ==========================================================================================
# Differences in microbial composition (beta-diversity) for all sample types across two GDD
# ==========================================================================================

# ================================
# read in the non-normalized data
# ================================

root <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/root_bac_sub_no_bad_Filtered_final_non_normalized.rds")
stalk <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/stalk_bac_sub_no_bad_Filtered_final_non_normalized.rds")
soil <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/soil_bac_sub_no_bad_Filtered_final_non_normalized.rds")
rhizo <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/rhizo_bac_sub_no_bad_Filtered_final_non_normalized.rds")


#for soil samples
ps_soil <- soil 
ps_tree_soil = phy_tree(ps_soil)

# Check if binary (dichotomy) & multifurcating (polytomy) trees
sprintf("Is tree binary: %s", is.binary(ps_tree_soil))

# If FASLE, randomly resolve polytomies and replace tree in "ps"
phy_tree(ps_soil) = multi2di(ps_tree_soil)
sprintf("Is tree binary: %s", is.binary(phy_tree(ps_soil)))

#reorder the facet grid 
ps_soil@sam_data[["Sample_Type"]] <- factor(ps_soil@sam_data[["Sample_Type"]], 
                                       levels=c("Soil", "Rhizosphere", "Root", "Stalk"))

# Convert phyloseq object to rbiom object
biom_soil <- as_rbiom(ps_soil)
biom_soil

#view th biom file
glimpse(biom_soil, width = NULL)

#calculate the unifrac distance between samples
biom_soil_dist <- bdiv_distmat(biom_soil, 'unifrac')
biom_soil_dist 

#run a Principal Coordinates Analysis on my distance matrix.
pcoa_biom_soil <- ape::pcoa(biom_soil_dist)
pcoa_biom_soil

#Calculate the Percent variance
PCoA1_biom_soil <- round(pcoa_biom_soil$values$Relative_eig[1] * 100, 2)
PCoA2_biom_soil <- round(pcoa_biom_soil$values$Relative_eig[2] * 100, 2)

#Ordinate samples based on beta diversity distances
wt.unifrac_soil <- bdiv_ord_plot(
  biom_soil,
  layers   = "pemt",
  bdiv     = "UniFrac",
  stat.by  = "Growing_Degree_Days",
  rank     = NULL
)

# Add points and title
wt.unifrac_soil <- wt.unifrac_soil +
  guides(shape = "none") +
  ggtitle("Weighted UniFrac Diversity PCoA - Soil") +
  geom_point(aes(color = Growing_Degree_Days), size = 6)

# Apply theme, legend, and color scale
wt.unifrac_soil <- wt.unifrac_soil +
  theme_gray() +
  theme(
    plot.title = element_text(family = "serif", face = "bold", size = 15),
    axis.title = element_text(family = "serif", size = 15, colour = "black", face = "bold"),
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 15, colour = "black"),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 15, colour = "black"),
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(size = 15, colour = "black"),
    legend.key.width = unit(1, "cm")) +
  xlab("PCoA 1 (12.32%)") +
  ylab("PCoA 2 (10.87%)") +
  scale_color_manual(values = colors_2, name = "Growing Degree Days")

# Draw plot
wt.unifrac_soil

#save plot as image
ggsave("wt.unifrac_soil.png", plot = wt.unifrac_soil , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save plot as pdf
ggsave("wt.unifrac_soil.pdf", plot = wt.unifrac_soil, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")



#for Rhizosphere samples
ps_rhizo <- rhizo 
ps_tree_rhizo = phy_tree(ps_rhizo)

# Check if binary (dichotomy) & multifurcating (polytomy) trees
sprintf("Is tree binary: %s", is.binary(ps_tree_rhizo))

# If FASLE, randomly resolve polytomies and replace tree in "ps"
phy_tree(ps_rhizo) = multi2di(ps_tree_rhizo)
sprintf("Is tree binary: %s", is.binary(phy_tree(ps_rhizo)))

#reorder the facet grid 
ps_rhizo@sam_data[["Sample_Type"]] <- factor(ps_rhizo@sam_data[["Sample_Type"]], 
                                            levels=c("Soil", "Rhizosphere", "Root", "Stalk"))

# Convert phyloseq object to rbiom object
biom_rhizo <- as_rbiom(ps_rhizo)
biom_rhizo

#view th biom file
glimpse(biom_rhizo, width = NULL)

#calculate the unifrac distance between samples
biom_rhizo_dist <- bdiv_distmat(biom_rhizo, 'unifrac')
biom_rhizo_dist 

#run a Principal Coordinates Analysis on my distance matrix.
pcoa_biom_rhizo <- ape::pcoa(biom_rhizo_dist)
pcoa_biom_rhizo

#Calculate the Percent variance
PCoA1_biom_rhizo <- round(pcoa_biom_rhizo$values$Relative_eig[1] * 100, 2)
PCoA2_biom_rhizo <- round(pcoa_biom_rhizo$values$Relative_eig[2] * 100, 2)

#Ordinate samples based on beta diversity distances
wt.unifrac_rhizo <- bdiv_ord_plot(
  biom_rhizo,
  layers   = "pemt",
  bdiv     = "UniFrac",
  stat.by  = "Growing_Degree_Days",
  rank     = NULL
)

# Add points and title
wt.unifrac_rhizo <- wt.unifrac_rhizo +
  guides(shape = "none") +
  ggtitle("Weighted UniFrac Diversity PCoA - Rhizosphere") +
  geom_point(aes(color = Growing_Degree_Days), size = 6)

# Apply theme, legend, and color scale
wt.unifrac_rhizo <- wt.unifrac_rhizo +
  theme_gray() +
  theme(
    plot.title = element_text(family = "serif", face = "bold", size = 15),
    axis.title = element_text(family = "serif", size = 15, colour = "black", face = "bold"),
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 15, colour = "black"),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 15, colour = "black"),
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(size = 15, colour = "black"),
    legend.key.width = unit(1, "cm")) +
  xlab("PCoA 1 (15.97%)") +
  ylab("PCoA 2 (10.14%)") +
  scale_color_manual(values = colors_2, name = "Growing Degree Days")

# Draw plot
wt.unifrac_rhizo


#save plot as image
ggsave("wt.unifrac_rhizo.png", plot = wt.unifrac_rhizo, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save plot as pdf
ggsave("wt.unifrac_rhizo.pdf", plot = wt.unifrac_rhizo, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")



#for Root samples
ps_root <- root 
ps_tree_root = phy_tree(ps_root)

# Check if binary (dichotomy) & multifurcating (polytomy) trees
sprintf("Is tree binary: %s", is.binary(ps_tree_root))

# If FASLE, randomly resolve polytomies and replace tree in "ps"
phy_tree(ps_root) = multi2di(ps_tree_root)
sprintf("Is tree binary: %s", is.binary(phy_tree(ps_root)))

#reorder the facet grid 
ps_root@sam_data[["Sample_Type"]] <- factor(ps_root@sam_data[["Sample_Type"]], 
                                             levels=c("Soil", "Rhizosphere", "Root", "Stalk"))

# Convert phyloseq object to rbiom object
biom_root <- as_rbiom(ps_root)
biom_root

#view th biom file
glimpse(biom_root, width = NULL)

#calculate the unifrac distance between samples
biom_root_dist <- bdiv_distmat(biom_root, 'unifrac')
biom_root_dist 

#run a Principal Coordinates Analysis on my distance matrix.
pcoa_biom_root <- ape::pcoa(biom_root_dist)
pcoa_biom_root

#Calculate the Percent variance
PCoA1_biom_root <- round(pcoa_biom_root$values$Relative_eig[1] * 100, 2)
PCoA2_biom_root <- round(pcoa_biom_root$values$Relative_eig[2] * 100, 2)

#Ordinate samples based on beta diversity distances
wt.unifrac_root <- bdiv_ord_plot(
  biom_root,
  layers   = "pemt",
  bdiv     = "UniFrac",
  stat.by  = "Growing_Degree_Days",
  rank     = NULL
)

# Add points and title
wt.unifrac_root <- wt.unifrac_root +
  guides(shape = "none") +
  ggtitle("Weighted UniFrac Diversity PCoA - Root") +
  geom_point(aes(color = Growing_Degree_Days), size = 6)

# Apply theme, legend, and color scale
wt.unifrac_root <- wt.unifrac_root +
  theme_gray() +
  theme(
    plot.title = element_text(family = "serif", face = "bold", size = 15),
    axis.title = element_text(family = "serif", size = 15, colour = "black", face = "bold"),
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 15, colour = "black"),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 15, colour = "black"),
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(size = 15, colour = "black"),
    legend.key.width = unit(1, "cm")) +
  xlab("PCoA 1 (20.95%)") +
  ylab("PCoA 2 (13.12%)") +
  scale_color_manual(values = colors_2, name = "Growing Degree Days")

# Draw plot
wt.unifrac_root


#save plot as image
ggsave("wt.unifrac_root.png", plot = wt.unifrac_root, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save plot as pdf
ggsave("wt.unifrac_root.pdf", plot = wt.unifrac_root, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")




#for Stalk samples
ps_stalk <- stalk 
ps_tree_stalk = phy_tree(ps_stalk)

# Check if binary (dichotomy) & multifurcating (polytomy) trees
sprintf("Is tree binary: %s", is.binary(ps_tree_stalk))

# If FASLE, randomly resolve polytomies and replace tree in "ps"
phy_tree(ps_stalk) = multi2di(ps_tree_stalk)
sprintf("Is tree binary: %s", is.binary(phy_tree(ps_stalk)))

#reorder the facet grid 
ps_stalk@sam_data[["Sample_Type"]] <- factor(ps_stalk@sam_data[["Sample_Type"]], 
                                            levels=c("Soil", "Rhizosphere", "Root", "Stalk"))

# Convert phyloseq object to rbiom object
biom_stalk <- as_rbiom(ps_stalk)
biom_stalk

#view th biom file
glimpse(biom_stalk, width = NULL)

#calculate the unifrac distance between samples
biom_stalk_dist <- bdiv_distmat(biom_stalk, 'unifrac')
biom_stalk_dist 

#run a Principal Coordinates Analysis on my distance matrix.
pcoa_biom_stalk <- ape::pcoa(biom_stalk_dist)
pcoa_biom_stalk

#Calculate the Percent variance
PCoA1_biom_stalk <- round(pcoa_biom_stalk$values$Relative_eig[1] * 100, 2)
PCoA2_biom_stalk <- round(pcoa_biom_stalk$values$Relative_eig[2] * 100, 2)

#Ordinate samples based on beta diversity distances
wt.unifrac_stalk <- bdiv_ord_plot(
  biom_stalk,
  layers   = "pemt",
  bdiv     = "UniFrac",
  stat.by  = "Growing_Degree_Days",
  rank     = NULL
)

# Add points and title
wt.unifrac_stalk <- wt.unifrac_stalk +
  guides(shape = "none") +
  ggtitle("Weighted UniFrac Diversity PCoA - Stalk") +
  geom_point(aes(color = Growing_Degree_Days), size = 6)

# Apply theme, legend, and color scale
wt.unifrac_stalk <- wt.unifrac_stalk +
  theme_gray() +
  theme(
    plot.title = element_text(family = "serif", face = "bold", size = 15),
    axis.title = element_text(family = "serif", size = 15, colour = "black", face = "bold"),
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 15, colour = "black"),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 15, colour = "black"),
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(size = 15, colour = "black"),
    legend.key.width = unit(1, "cm")) +
  xlab("PCoA 1 (26.77%)") +
  ylab("PCoA 2 (13.14%)") +
  scale_color_manual(values = colors_2, name = "Growing Degree Days")

# Draw plot
wt.unifrac_stalk


#save plot as image
ggsave("wt.unifrac_stalk.png", plot = wt.unifrac_stalk, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save plot as pdf
ggsave("wt.unifrac_stalk.pdf", plot = wt.unifrac_stalk, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")

