#Load required libraries
library(phyloseq)
library(geosphere)
library(ggplot2)
library(ggpubr)
library(ape)
library(rbiom)

set.seed(1234)

# ================================
# read in the non-normalized data
# ================================
root_bac <- readRDS("C:/Aduragbemi/Microbiome/RDS/root_bac_sub_no_bad_Filtered_final_non_normalized.rds")
stalk_bac <- readRDS("C:/Aduragbemi/Microbiome/RDS/stalk_bac_sub_no_bad_Filtered_final_non_normalized.rds")
soil_bac  <- readRDS("C:/Aduragbemi/Microbiome/RDS/soil_bac_sub_no_bad_Filtered_final_non_normalized.rds")
rhizo_bac <- readRDS("C:/Aduragbemi/Microbiome/RDS/rhizo_bac_sub_no_bad_Filtered_final_non_normalized.rds")


# ================================================================================
# correlation between geographic distance  (x axis) and weighted Unifrac distance
# ================================================================================

#For Soil Sample
# Assign to shorter variable name
ps_soil <- soil_bac
ps_soil_tree = phy_tree(ps_soil)

# Check if binary (dichotomy) & multifurcating (polytomy) trees
sprintf("Is tree binary: %s", is.binary(ps_soil_tree))

# If FASLE, randomly resolve polytomies and replace tree in "ps"
phy_tree(ps_soil) = multi2di(ps_soil_tree)
sprintf("Is tree binary: %s", is.binary(phy_tree(ps_soil)))

# Convert phyloseq object to rbiom object
biom_soil <- as_rbiom(ps_soil)
biom_soil

#view th biom file
glimpse(biom_soil, width = NULL)

#calculate the unifrac distance between samples
biom_dm_soil <- bdiv_distmat(biom_soil, 'unifrac')
biom_dm_soil

#Extract GPS Coordinates
coords_soil <- sample_data(ps_soil)[, c("GPSlatitude", "GPSlongitude")]
coords_soil <- as.data.frame(coords_soil)
coords_soil <- as.matrix(coords_soil)

#Compute Physical Distances (Haversine in meters) 
phys_dist_matrix_soil <- distm(coords_soil, fun = distHaversine)
phys_dist_soil <- as.dist(phys_dist_matrix_soil)

#Convert to Vectors for Correlation
unifrac_vec_soil <- as.vector(as.dist(biom_dm_soil))
phys_vec_soil <- as.vector(phys_dist_soil)

#Log-Transform Physical Distances
phys_vec_log_soil <- log1p(phys_vec_soil)  # log(1 + distance)

#Spearman Correlation
cor_result_soil <- cor.test(unifrac_vec_soil, phys_vec_log_soil, method = "spearman")
print(cor_result_soil)

#Create Data Frame for Plotting
df_dist_soil <- data.frame(
  UniFrac_Distance_soil = unifrac_vec_soil,
  Physical_Distance_log_soil = phys_vec_log_soil
)

# Apply theme, legend, and color scale and draw the plot
plot_soil <- ggplot(df_dist_soil, aes(x = Physical_Distance_log_soil, y = UniFrac_Distance_soil)) +
  theme_gray(base_size = 14) +  # Gray theme applied here
  geom_point(alpha = 0.4, size = 6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
  stat_cor(method = "spearman", 
           label.x.npc = "left", label.y.npc = "top", 
           size = 15, 
           label.sep = ", ") +
  labs(
    title = "Soil microbiome - Spearman Correlation",
    x = "Geographic Distance (m)",
    y = "Weighted UniFrac Distance"
  ) +
  theme(
    plot.title = element_text(family = "serif", face = "bold", size = 15),
    axis.title = element_text(family = "serif", size = 15, colour = "black", face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 15, colour = "black"),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 15, colour = "black"),
    legend.position = "bottom",
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(size = 15, colour = "black")
  )

plot_soil

#save plot as image
ggsave("plot_soil.png", plot = plot_soil , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save plot as pdf
ggsave("plot_soil.pdf", plot = plot_soil, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")




#For Rhizosphere Sample
# Assign to shorter variable name
ps_rhizo <- rhizo_bac
ps_rhizo_tree = phy_tree(ps_rhizo)

# Check if binary (dichotomy) & multifurcating (polytomy) trees
sprintf("Is tree binary: %s", is.binary(ps_rhizo_tree))

# If FASLE, randomly resolve polytomies and replace tree in "ps"
phy_tree(ps_rhizo) = multi2di(ps_rhizo_tree)
sprintf("Is tree binary: %s", is.binary(phy_tree(ps_rhizo)))

# Convert phyloseq object to rbiom object
biom_rhizo <- as_rbiom(ps_rhizo)
biom_rhizo

#view th biom file
glimpse(biom_rhizo, width = NULL)

#calculate the unifrac distance between samples
biom_dm_rhizo <- bdiv_distmat(biom_rhizo, 'unifrac')
biom_dm_rhizo

#Extract GPS Coordinates
coords_rhizo <- sample_data(ps_rhizo)[, c("GPSlatitude", "GPSlongitude")]
coords_rhizo <- as.data.frame(coords_rhizo)
coords_rhizo <- as.matrix(coords_rhizo)

#Compute Physical Distances (Haversine in meters) 
phys_dist_matrix_rhizo <- distm(coords_rhizo, fun = distHaversine)
phys_dist_rhizo <- as.dist(phys_dist_matrix_rhizo)

#Convert to Vectors for Correlation
unifrac_vec_rhizo <- as.vector(as.dist(biom_dm_rhizo))
phys_vec_rhizo <- as.vector(phys_dist_rhizo)

#Log-Transform Physical Distances
phys_vec_log_rhizo <- log1p(phys_vec_rhizo)  # log(1 + distance)

#Spearman Correlation
cor_result_rhizo <- cor.test(unifrac_vec_rhizo, phys_vec_log_rhizo, method = "spearman")
print(cor_result_rhizo)

#Create Data Frame for Plotting
df_dist_rhizo <- data.frame(
  UniFrac_Distance_rhizo = unifrac_vec_rhizo,
  Physical_Distance_log_rhizo = phys_vec_log_rhizo
)

# Apply theme, legend, and color scale and draw the plot
plot_rhizo <- ggplot(df_dist_rhizo, aes(x = Physical_Distance_log_rhizo, y = UniFrac_Distance_rhizo)) +
  theme_gray(base_size = 14) +  # Gray theme applied here
  geom_point(alpha = 0.4, size = 6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
  stat_cor(method = "spearman", 
           label.x.npc = "left", label.y.npc = "top", 
           size = 15, 
           label.sep = ", ") +
  labs(
    title = "Rhizosphere microbiome - Spearman Correlation",
    x = "Geographic Distance (m)",
    y = "Weighted UniFrac Distance"
  ) +
  theme(
    plot.title = element_text(family = "serif", face = "bold", size = 15),
    axis.title = element_text(family = "serif", size = 15, colour = "black", face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 15, colour = "black"),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 15, colour = "black"),
    legend.position = "bottom",
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(size = 15, colour = "black")
  )

plot_rhizo

#save plot as image
ggsave("plot_rhizo.png", plot = plot_rhizo , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save plot as pdf
ggsave("plot_rhizo.pdf", plot = plot_rhizo, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")




#For Root Sample
# Assign to shorter variable name
ps_root <- root_bac
ps_root_tree = phy_tree(ps_root)

# Check if binary (dichotomy) & multifurcating (polytomy) trees
sprintf("Is tree binary: %s", is.binary(ps_root_tree))

# If FASLE, randomly resolve polytomies and replace tree in "ps"
phy_tree(ps_root) = multi2di(ps_root_tree)
sprintf("Is tree binary: %s", is.binary(phy_tree(ps_root)))

# Convert phyloseq object to rbiom object
biom_root <- as_rbiom(ps_root)
biom_root

#view th biom file
glimpse(biom_root, width = NULL)

#calculate the unifrac distance between samples
biom_dm_root <- bdiv_distmat(biom_root, 'unifrac')
biom_dm_root

#Extract GPS Coordinates
coords_root <- sample_data(ps_root)[, c("GPSlatitude", "GPSlongitude")]
coords_root <- as.data.frame(coords_root)
coords_root <- as.matrix(coords_root)

#Compute Physical Distances (Haversine in meters) 
phys_dist_matrix_root <- distm(coords_root, fun = distHaversine)
phys_dist_root <- as.dist(phys_dist_matrix_root)

#Convert to Vectors for Correlation
unifrac_vec_root <- as.vector(as.dist(biom_dm_root))
phys_vec_root <- as.vector(phys_dist_root)

#Log-Transform Physical Distances
phys_vec_log_root <- log1p(phys_vec_root)  # log(1 + distance)

#Spearman Correlation
cor_result_root <- cor.test(unifrac_vec_root, phys_vec_log_root, method = "spearman")
print(cor_result_root)

#Create Data Frame for Plotting
df_dist_root <- data.frame(
  UniFrac_Distance_root = unifrac_vec_root,
  Physical_Distance_log_root = phys_vec_log_root
)

# Apply theme, legend, and color scale and draw the plot
plot_root <- ggplot(df_dist_root, aes(x = Physical_Distance_log_root, y = UniFrac_Distance_root)) +
  theme_gray(base_size = 14) +  # Gray theme applied here
  geom_point(alpha = 0.4, size = 6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
  stat_cor(method = "spearman", 
           label.x.npc = "left", label.y.npc = "top", 
           size = 15, 
           label.sep = ", ") +
  labs(
    title = "Root microbiome - Spearman Correlation",
    x = "Geographic Distance (m)",
    y = "Weighted UniFrac Distance"
  ) +
  theme(
    plot.title = element_text(family = "serif", face = "bold", size = 15),
    axis.title = element_text(family = "serif", size = 15, colour = "black", face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 15, colour = "black"),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 15, colour = "black"),
    legend.position = "bottom",
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(size = 15, colour = "black")
  )

plot_root

#save plot as image
ggsave("plot_root.png", plot = plot_root , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save plot as pdf
ggsave("plot_root.pdf", plot = plot_root, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")







#For Stalk Sample
#Assign to shorter variable name
ps_stalk <- stalk_bac
ps_stalk_tree = phy_tree(ps_stalk)

# Check if binary (dichotomy) & multifurcating (polytomy) trees
sprintf("Is tree binary: %s", is.binary(ps_stalk_tree))

# If FASLE, randomly resolve polytomies and replace tree in "ps"
phy_tree(ps_stalk) = multi2di(ps_stalk_tree)
sprintf("Is tree binary: %s", is.binary(phy_tree(ps_stalk)))

# Convert phyloseq object to rbiom object
biom_stalk <- as_rbiom(ps_stalk)
biom_stalk

#view th biom file
glimpse(biom_stalk, width = NULL)

#calculate the unifrac distance between samples
biom_dm_stalk <- bdiv_distmat(biom_stalk, 'unifrac')
biom_dm_stalk

#Extract GPS Coordinates
coords_stalk <- sample_data(ps_stalk)[, c("GPSlatitude", "GPSlongitude")]
coords_stalk <- as.data.frame(coords_stalk)
coords_stalk <- as.matrix(coords_stalk)

#Compute Physical Distances (Haversine in meters) 
phys_dist_matrix_stalk <- distm(coords_stalk, fun = distHaversine)
phys_dist_stalk <- as.dist(phys_dist_matrix_stalk)

#Convert to Vectors for Correlation
unifrac_vec_stalk <- as.vector(as.dist(biom_dm_stalk))
phys_vec_stalk <- as.vector(phys_dist_stalk)

#Log-Transform Physical Distances
phys_vec_log_stalk <- log1p(phys_vec_stalk)  # log(1 + distance)

#Spearman Correlation
cor_result_stalk <- cor.test(unifrac_vec_stalk, phys_vec_log_stalk, method = "spearman")
print(cor_result_stalk)

#Create Data Frame for Plotting
df_dist_stalk <- data.frame(
  UniFrac_Distance_stalk = unifrac_vec_stalk,
  Physical_Distance_log_stalk = phys_vec_log_stalk
)

# Apply theme, legend, and color scale and draw the plot
plot_stalk <- ggplot(df_dist_stalk, aes(x = Physical_Distance_log_stalk, y = UniFrac_Distance_stalk)) +
  theme_gray(base_size = 14) +  # Gray theme applied here
  geom_point(alpha = 0.4, size = 6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
  stat_cor(method = "spearman", 
           label.x.npc = "left", label.y.npc = "top", 
           size = 15, 
           label.sep = ", ") +
  labs(
    title = "Stalk microbiome - Spearman Correlation",
    x = "Geographic Distance (m)",
    y = "Weighted UniFrac Distance"
  ) +
  theme(
    plot.title = element_text(family = "serif", face = "bold", size = 15),
    axis.title = element_text(family = "serif", size = 15, colour = "black", face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 15, colour = "black"),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 15, colour = "black"),
    legend.position = "bottom",
    legend.title = element_text(colour = "black", size = 15, face = "bold"),
    legend.text = element_text(size = 15, colour = "black")
  )

plot_stalk

#save plot as image
ggsave("plot_stalk.png", plot = plot_stalk , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save plot as pdf
ggsave("plot_stalk.pdf", plot = plot_stalk, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")

