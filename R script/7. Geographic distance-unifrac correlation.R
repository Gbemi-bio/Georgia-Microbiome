library(phyloseq)
library(geosphere)
library(ggplot2)
library(ggpubr)
library(ape)
library(rbiom)
library(vegan)


set.seed(1234)

#===============================
# Read data
#===============================

soil_bac <- readRDS("C:/Aduragbemi/Microbiome/RDS/soil_bac_sub_no_bad_Filtered_final_non_normalized.rds")
rhizo_bac <- readRDS("C:/Aduragbemi/Microbiome/RDS/rhizo_bac_sub_no_bad_Filtered_final_non_normalized.rds")
root_bac <- readRDS("C:/Aduragbemi/Microbiome/RDS/root_bac_sub_no_bad_Filtered_final_non_normalized.rds")
stalk_bac <- readRDS("C:/Aduragbemi/Microbiome/RDS/stalk_bac_sub_no_bad_Filtered_final_non_normalized.rds")


#For soil samples

#==========================================================
# Soil samples - 600 GDD
#==========================================================

ps_soil_600 <- subset_samples(
  soil_bac,
  Growing_Degree_Days == "600"
)

soil_tree_600 <- phy_tree(ps_soil_600)

# Check if tree is binary
sprintf("Is tree binary: %s", is.binary(soil_tree_600))

# Resolve polytomies
phy_tree(ps_soil_600) <- multi2di(soil_tree_600)

sprintf("Is tree binary: %s", is.binary(phy_tree(ps_soil_600)))


# Convert phyloseq object to rbiom object
biom_soil_600 <- as_rbiom(ps_soil_600)
biom_soil_600

# View the biom file
glimpse(biom_soil_600, width = NULL)


# Calculate UniFrac distance between samples
biom_dm_soil_600 <- bdiv_distmat(
  biom_soil_600,
  bdiv = "UniFrac",
  weighted = TRUE,
  normalized = TRUE
)
biom_dm_soil_600


# Convert UniFrac distance to matrix
biom_dm_soil_600 <- as.matrix(biom_dm_soil_600)


#-------------------------------------------------
# Check sample IDs
#-------------------------------------------------

# Check if samples are the same
setequal(
  sample_names(ps_soil_600),
  rownames(biom_dm_soil_600)
)

# Check if sample order is the same
identical(
  sample_names(ps_soil_600),
  rownames(biom_dm_soil_600)
)

# Check UniFrac row and column names
identical(
  rownames(biom_dm_soil_600),
  colnames(biom_dm_soil_600)
)


#-------------------------------------------------
# Extract GPS coordinates
#-------------------------------------------------

# Longitude first, latitude second
coords_soil_600 <- sample_data(ps_soil_600)[,
                                            c("GPSlongitude", "GPSlatitude")
]

coords_soil_600 <- as.data.frame(coords_soil_600)

# Convert coordinates to numeric
coords_soil_600$GPSlongitude <- as.numeric(
  as.character(coords_soil_600$GPSlongitude)
)

coords_soil_600$GPSlatitude <- as.numeric(
  as.character(coords_soil_600$GPSlatitude)
)


# Check whether sample IDs match
setequal(
  rownames(coords_soil_600),
  rownames(biom_dm_soil_600)
)


# Match coordinate order with UniFrac matrix
coords_soil_600 <- coords_soil_600[
  rownames(biom_dm_soil_600),
  ,
  drop = FALSE
]

# Match UniFrac column order with row order
biom_dm_soil_600 <- biom_dm_soil_600[
  rownames(biom_dm_soil_600),
  rownames(biom_dm_soil_600)
]


# Verify sample order
identical(
  rownames(coords_soil_600),
  rownames(biom_dm_soil_600)
)


# Check GPS coordinates
summary(coords_soil_600)

colSums(is.na(coords_soil_600))


# Convert coordinates to matrix
coords_soil_600 <- as.matrix(coords_soil_600)


#-------------------------------------------------
# Compute Physical Distances (Haversine)
#-------------------------------------------------

# Calculate geographic distance in meters
phys_dist_matrix_soil_600 <- distm(
  coords_soil_600,
  fun = distHaversine
)

# Convert meters to kilometers
phys_dist_matrix_soil_600 <- phys_dist_matrix_soil_600 / 1000

phys_dist_soil_600 <- as.dist(phys_dist_matrix_soil_600)


#-------------------------------------------------
# Convert to Vectors for Correlation
#-------------------------------------------------

unifrac_vec_soil_600 <- as.vector(
  as.dist(biom_dm_soil_600)
)

phys_vec_soil_600 <- as.vector(
  phys_dist_soil_600
)


#-------------------------------------------------
# Check Geographic Distances (km)
#-------------------------------------------------

summary(phys_vec_soil_600)

max(phys_vec_soil_600)


#-------------------------------------------------
# Log-Transform Physical Distances
#-------------------------------------------------

phys_vec_log_soil_600 <- log1p(phys_vec_soil_600)


#-------------------------------------------------
# Spearman Correlation
#-------------------------------------------------

cor_result_soil_600 <- cor.test(
  unifrac_vec_soil_600,
  phys_vec_log_soil_600,
  method = "spearman",
  exact = FALSE
)

print(cor_result_soil_600)


#-------------------------------------------------
# Create Data Frame for Plotting
#-------------------------------------------------

df_dist_soil_600 <- data.frame(
  UniFrac_Distance_soil_600 = unifrac_vec_soil_600,
  Physical_Distance_km_soil_600 = phys_vec_soil_600,
  Physical_Distance_log_soil_600 = phys_vec_log_soil_600
)

head(df_dist_soil_600)


#------------------------------------------------------------------------------------------------------------
# Perform the Mantel test, Since you are analyzing geographic distance against weighted UniFrac distance
#---------------------------------------------------------------------------------------------------------
mantel_soil_600 <- mantel(
  biom_dm_soil_600,
  as.matrix(phys_dist_soil_600),
  method = "spearman",
  permutations = 999
)

print(mantel_soil_600)



#==========================================================
# Soil samples - 1400 GDD
#==========================================================

ps_soil_1400 <- subset_samples(
  soil_bac,
  Growing_Degree_Days == "1400"
)

soil_tree_1400 <- phy_tree(ps_soil_1400)

# Check if tree is binary
sprintf("Is tree binary: %s", is.binary(soil_tree_1400))

# Resolve polytomies
phy_tree(ps_soil_1400) <- multi2di(soil_tree_1400)

sprintf("Is tree binary: %s", is.binary(phy_tree(ps_soil_1400)))


# Convert phyloseq object to rbiom object
biom_soil_1400 <- as_rbiom(ps_soil_1400)
biom_soil_1400

# View the biom file
glimpse(biom_soil_1400, width = NULL)


# Calculate UniFrac distance between samples
biom_dm_soil_1400 <- bdiv_distmat(
  biom_soil_1400,
  bdiv = "UniFrac",
  weighted = TRUE,
  normalized = TRUE
)
biom_dm_soil_1400


# Convert UniFrac distance to matrix
biom_dm_soil_1400 <- as.matrix(biom_dm_soil_1400)


#-------------------------------------------------
# Check sample IDs
#-------------------------------------------------

# Check if samples are the same
setequal(
  sample_names(ps_soil_1400),
  rownames(biom_dm_soil_1400)
)

# Check if sample order is the same
identical(
  sample_names(ps_soil_1400),
  rownames(biom_dm_soil_1400)
)

# Check UniFrac row and column names
identical(
  rownames(biom_dm_soil_1400),
  colnames(biom_dm_soil_1400)
)


#-------------------------------------------------
# Extract GPS coordinates
#-------------------------------------------------

# Longitude first, latitude second
coords_soil_1400 <- sample_data(ps_soil_1400)[,
                                            c("GPSlongitude", "GPSlatitude")
]

coords_soil_1400 <- as.data.frame(coords_soil_1400)

# Convert coordinates to numeric
coords_soil_1400$GPSlongitude <- as.numeric(
  as.character(coords_soil_1400$GPSlongitude)
)

coords_soil_1400$GPSlatitude <- as.numeric(
  as.character(coords_soil_1400$GPSlatitude)
)


# Check whether sample IDs match
setequal(
  rownames(coords_soil_1400),
  rownames(biom_dm_soil_1400)
)


# Match coordinate order with UniFrac matrix
coords_soil_1400 <- coords_soil_1400[
  rownames(biom_dm_soil_1400),
  ,
  drop = FALSE
]

# Match UniFrac column order with row order
biom_dm_soil_1400 <- biom_dm_soil_1400[
  rownames(biom_dm_soil_1400),
  rownames(biom_dm_soil_1400)
]


# Verify sample order
identical(
  rownames(coords_soil_1400),
  rownames(biom_dm_soil_1400)
)


# Check GPS coordinates
summary(coords_soil_1400)

colSums(is.na(coords_soil_1400))


# Convert coordinates to matrix
coords_soil_1400 <- as.matrix(coords_soil_1400)


#-------------------------------------------------
# Compute Physical Distances (Haversine)
#-------------------------------------------------

# Calculate geographic distance in meters
phys_dist_matrix_soil_1400 <- distm(
  coords_soil_1400,
  fun = distHaversine
)

# Convert meters to kilometers
phys_dist_matrix_soil_1400 <- phys_dist_matrix_soil_1400 / 1000

phys_dist_soil_1400 <- as.dist(phys_dist_matrix_soil_1400)


#-------------------------------------------------
# Convert to Vectors for Correlation
#-------------------------------------------------

unifrac_vec_soil_1400 <- as.vector(
  as.dist(biom_dm_soil_1400)
)

phys_vec_soil_1400 <- as.vector(
  phys_dist_soil_1400
)


#-------------------------------------------------
# Check Geographic Distances (km)
#-------------------------------------------------

summary(phys_vec_soil_1400)

max(phys_vec_soil_1400)


#-------------------------------------------------
# Log-Transform Physical Distances
#-------------------------------------------------

phys_vec_log_soil_1400 <- log1p(phys_vec_soil_1400)


#-------------------------------------------------
# Spearman Correlation
#-------------------------------------------------

cor_result_soil_1400 <- cor.test(
  unifrac_vec_soil_1400,
  phys_vec_log_soil_1400,
  method = "spearman",
  exact = FALSE
)

print(cor_result_soil_1400)


#-------------------------------------------------
# Create Data Frame for Plotting
#-------------------------------------------------

df_dist_soil_1400 <- data.frame(
  UniFrac_Distance_soil_1400 = unifrac_vec_soil_1400,
  Physical_Distance_km_soil_1400 = phys_vec_soil_1400,
  Physical_Distance_log_soil_1400 = phys_vec_log_soil_1400
)

head(df_dist_soil_1400)


#------------------------------------------------------------------------------------------------------------
# Perform the Mantel test, Since you are analyzing geographic distance against weighted UniFrac distance
#---------------------------------------------------------------------------------------------------------
mantel_soil_1400 <- mantel(
  biom_dm_soil_1400,
  as.matrix(phys_dist_soil_1400),
  method = "spearman",
  permutations = 999
)

print(mantel_soil_1400)





#For Rhizosphere samples


#==========================================================
# rhizo samples - 600 GDD
#==========================================================

ps_rhizo_600 <- subset_samples(
  rhizo_bac,
  Growing_Degree_Days == "600"
)

rhizo_tree_600 <- phy_tree(ps_rhizo_600)

# Check if tree is binary
sprintf("Is tree binary: %s", is.binary(rhizo_tree_600))

# Resolve polytomies
phy_tree(ps_rhizo_600) <- multi2di(rhizo_tree_600)

sprintf("Is tree binary: %s", is.binary(phy_tree(ps_rhizo_600)))


# Convert phyloseq object to rbiom object
biom_rhizo_600 <- as_rbiom(ps_rhizo_600)
biom_rhizo_600

# View the biom file
glimpse(biom_rhizo_600, width = NULL)


# Calculate UniFrac distance between samples
biom_dm_rhizo_600 <- bdiv_distmat(
  biom_rhizo_600,
  bdiv = "UniFrac",
  weighted = TRUE,
  normalized = TRUE
)
biom_dm_rhizo_600


# Convert UniFrac distance to matrix
biom_dm_rhizo_600 <- as.matrix(biom_dm_rhizo_600)


#-------------------------------------------------
# Check sample IDs
#-------------------------------------------------

# Check if samples are the same
setequal(
  sample_names(ps_rhizo_600),
  rownames(biom_dm_rhizo_600)
)

# Check if sample order is the same
identical(
  sample_names(ps_rhizo_600),
  rownames(biom_dm_rhizo_600)
)

# Check UniFrac row and column names
identical(
  rownames(biom_dm_rhizo_600),
  colnames(biom_dm_rhizo_600)
)


#-------------------------------------------------
# Extract GPS coordinates
#-------------------------------------------------

# Longitude first, latitude second
coords_rhizo_600 <- sample_data(ps_rhizo_600)[,
                                            c("GPSlongitude", "GPSlatitude")
]

coords_rhizo_600 <- as.data.frame(coords_rhizo_600)

# Convert coordinates to numeric
coords_rhizo_600$GPSlongitude <- as.numeric(
  as.character(coords_rhizo_600$GPSlongitude)
)

coords_rhizo_600$GPSlatitude <- as.numeric(
  as.character(coords_rhizo_600$GPSlatitude)
)


# Check whether sample IDs match
setequal(
  rownames(coords_rhizo_600),
  rownames(biom_dm_rhizo_600)
)


# Match coordinate order with UniFrac matrix
coords_rhizo_600 <- coords_rhizo_600[
  rownames(biom_dm_rhizo_600),
  ,
  drop = FALSE
]

# Match UniFrac column order with row order
biom_dm_rhizo_600 <- biom_dm_rhizo_600[
  rownames(biom_dm_rhizo_600),
  rownames(biom_dm_rhizo_600)
]


# Verify sample order
identical(
  rownames(coords_rhizo_600),
  rownames(biom_dm_rhizo_600)
)


# Check GPS coordinates
summary(coords_rhizo_600)

colSums(is.na(coords_rhizo_600))


# Convert coordinates to matrix
coords_rhizo_600 <- as.matrix(coords_rhizo_600)


#-------------------------------------------------
# Compute Physical Distances (Haversine)
#-------------------------------------------------

# Calculate geographic distance in meters
phys_dist_matrix_rhizo_600 <- distm(
  coords_rhizo_600,
  fun = distHaversine
)

# Convert meters to kilometers
phys_dist_matrix_rhizo_600 <- phys_dist_matrix_rhizo_600 / 1000

phys_dist_rhizo_600 <- as.dist(phys_dist_matrix_rhizo_600)


#-------------------------------------------------
# Convert to Vectors for Correlation
#-------------------------------------------------

unifrac_vec_rhizo_600 <- as.vector(
  as.dist(biom_dm_rhizo_600)
)

phys_vec_rhizo_600 <- as.vector(
  phys_dist_rhizo_600
)


#-------------------------------------------------
# Check Geographic Distances (km)
#-------------------------------------------------

summary(phys_vec_rhizo_600)

max(phys_vec_rhizo_600)


#-------------------------------------------------
# Log-Transform Physical Distances
#-------------------------------------------------

phys_vec_log_rhizo_600 <- log1p(phys_vec_rhizo_600)


#-------------------------------------------------
# Spearman Correlation
#-------------------------------------------------

cor_result_rhizo_600 <- cor.test(
  unifrac_vec_rhizo_600,
  phys_vec_log_rhizo_600,
  method = "spearman",
  exact = FALSE
)

print(cor_result_rhizo_600)


#-------------------------------------------------
# Create Data Frame for Plotting
#-------------------------------------------------

df_dist_rhizo_600 <- data.frame(
  UniFrac_Distance_rhizo_600 = unifrac_vec_rhizo_600,
  Physical_Distance_km_rhizo_600 = phys_vec_rhizo_600,
  Physical_Distance_log_rhizo_600 = phys_vec_log_rhizo_600
)

head(df_dist_rhizo_600)


#------------------------------------------------------------------------------------------------------------
# Perform the Mantel test, Since you are analyzing geographic distance against weighted UniFrac distance
#---------------------------------------------------------------------------------------------------------
mantel_rhizo_600 <- mantel(
  biom_dm_rhizo_600,
  as.matrix(phys_dist_rhizo_600),
  method = "spearman",
  permutations = 999
)

print(mantel_rhizo_600)



#==========================================================
# rhizo samples - 1400 GDD
#==========================================================

ps_rhizo_1400 <- subset_samples(
  rhizo_bac,
  Growing_Degree_Days == "1400"
)

rhizo_tree_1400 <- phy_tree(ps_rhizo_1400)

# Check if tree is binary
sprintf("Is tree binary: %s", is.binary(rhizo_tree_1400))

# Resolve polytomies
phy_tree(ps_rhizo_1400) <- multi2di(rhizo_tree_1400)

sprintf("Is tree binary: %s", is.binary(phy_tree(ps_rhizo_1400)))


# Convert phyloseq object to rbiom object
biom_rhizo_1400 <- as_rbiom(ps_rhizo_1400)
biom_rhizo_1400

# View the biom file
glimpse(biom_rhizo_1400, width = NULL)


# Calculate UniFrac distance between samples
biom_dm_rhizo_1400 <- bdiv_distmat(
  biom_rhizo_1400,
  bdiv = "UniFrac",
  weighted = TRUE,
  normalized = TRUE
)
biom_dm_rhizo_1400


# Convert UniFrac distance to matrix
biom_dm_rhizo_1400 <- as.matrix(biom_dm_rhizo_1400)


#-------------------------------------------------
# Check sample IDs
#-------------------------------------------------

# Check if samples are the same
setequal(
  sample_names(ps_rhizo_1400),
  rownames(biom_dm_rhizo_1400)
)

# Check if sample order is the same
identical(
  sample_names(ps_rhizo_1400),
  rownames(biom_dm_rhizo_1400)
)

# Check UniFrac row and column names
identical(
  rownames(biom_dm_rhizo_1400),
  colnames(biom_dm_rhizo_1400)
)


#-------------------------------------------------
# Extract GPS coordinates
#-------------------------------------------------

# Longitude first, latitude second
coords_rhizo_1400 <- sample_data(ps_rhizo_1400)[,
                                              c("GPSlongitude", "GPSlatitude")
]

coords_rhizo_1400 <- as.data.frame(coords_rhizo_1400)

# Convert coordinates to numeric
coords_rhizo_1400$GPSlongitude <- as.numeric(
  as.character(coords_rhizo_1400$GPSlongitude)
)

coords_rhizo_1400$GPSlatitude <- as.numeric(
  as.character(coords_rhizo_1400$GPSlatitude)
)


# Check whether sample IDs match
setequal(
  rownames(coords_rhizo_1400),
  rownames(biom_dm_rhizo_1400)
)


# Match coordinate order with UniFrac matrix
coords_rhizo_1400 <- coords_rhizo_1400[
  rownames(biom_dm_rhizo_1400),
  ,
  drop = FALSE
]

# Match UniFrac column order with row order
biom_dm_rhizo_1400 <- biom_dm_rhizo_1400[
  rownames(biom_dm_rhizo_1400),
  rownames(biom_dm_rhizo_1400)
]


# Verify sample order
identical(
  rownames(coords_rhizo_1400),
  rownames(biom_dm_rhizo_1400)
)


# Check GPS coordinates
summary(coords_rhizo_1400)

colSums(is.na(coords_rhizo_1400))


# Convert coordinates to matrix
coords_rhizo_1400 <- as.matrix(coords_rhizo_1400)


#-------------------------------------------------
# Compute Physical Distances (Haversine)
#-------------------------------------------------

# Calculate geographic distance in meters
phys_dist_matrix_rhizo_1400 <- distm(
  coords_rhizo_1400,
  fun = distHaversine
)

# Convert meters to kilometers
phys_dist_matrix_rhizo_1400 <- phys_dist_matrix_rhizo_1400 / 1000

phys_dist_rhizo_1400 <- as.dist(phys_dist_matrix_rhizo_1400)


#-------------------------------------------------
# Convert to Vectors for Correlation
#-------------------------------------------------

unifrac_vec_rhizo_1400 <- as.vector(
  as.dist(biom_dm_rhizo_1400)
)

phys_vec_rhizo_1400 <- as.vector(
  phys_dist_rhizo_1400
)


#-------------------------------------------------
# Check Geographic Distances (km)
#-------------------------------------------------

summary(phys_vec_rhizo_1400)

max(phys_vec_rhizo_1400)


#-------------------------------------------------
# Log-Transform Physical Distances
#-------------------------------------------------

phys_vec_log_rhizo_1400 <- log1p(phys_vec_rhizo_1400)


#-------------------------------------------------
# Spearman Correlation
#-------------------------------------------------

cor_result_rhizo_1400 <- cor.test(
  unifrac_vec_rhizo_1400,
  phys_vec_log_rhizo_1400,
  method = "spearman",
  exact = FALSE
)

print(cor_result_rhizo_1400)


#-------------------------------------------------
# Create Data Frame for Plotting
#-------------------------------------------------

df_dist_rhizo_1400 <- data.frame(
  UniFrac_Distance_rhizo_1400 = unifrac_vec_rhizo_1400,
  Physical_Distance_km_rhizo_1400 = phys_vec_rhizo_1400,
  Physical_Distance_log_rhizo_1400 = phys_vec_log_rhizo_1400
)

head(df_dist_rhizo_1400)


#------------------------------------------------------------------------------------------------------------
# Perform the Mantel test, Since you are analyzing geographic distance against weighted UniFrac distance
#---------------------------------------------------------------------------------------------------------
mantel_rhizo_1400 <- mantel(
  biom_dm_rhizo_1400,
  as.matrix(phys_dist_rhizo_1400),
  method = "spearman",
  permutations = 999
)

print(mantel_rhizo_1400)



#For Root samples

#==========================================================
# root samples - 600 GDD
#==========================================================

ps_root_600 <- subset_samples(
  root_bac,
  Growing_Degree_Days == "600"
)

root_tree_600 <- phy_tree(ps_root_600)

# Check if tree is binary
sprintf("Is tree binary: %s", is.binary(root_tree_600))

# Resolve polytomies
phy_tree(ps_root_600) <- multi2di(root_tree_600)

sprintf("Is tree binary: %s", is.binary(phy_tree(ps_root_600)))


# Convert phyloseq object to rbiom object
biom_root_600 <- as_rbiom(ps_root_600)
biom_root_600

# View the biom file
glimpse(biom_root_600, width = NULL)


# Calculate UniFrac distance between samples
biom_dm_root_600 <- bdiv_distmat(
  biom_root_600,
  bdiv = "UniFrac",
  weighted = TRUE,
  normalized = TRUE
)
biom_dm_root_600


# Convert UniFrac distance to matrix
biom_dm_root_600 <- as.matrix(biom_dm_root_600)


#-------------------------------------------------
# Check sample IDs
#-------------------------------------------------

# Check if samples are the same
setequal(
  sample_names(ps_root_600),
  rownames(biom_dm_root_600)
)

# Check if sample order is the same
identical(
  sample_names(ps_root_600),
  rownames(biom_dm_root_600)
)

# Check UniFrac row and column names
identical(
  rownames(biom_dm_root_600),
  colnames(biom_dm_root_600)
)


#-------------------------------------------------
# Extract GPS coordinates
#-------------------------------------------------

# Longitude first, latitude second
coords_root_600 <- sample_data(ps_root_600)[,
                                              c("GPSlongitude", "GPSlatitude")
]

coords_root_600 <- as.data.frame(coords_root_600)

# Convert coordinates to numeric
coords_root_600$GPSlongitude <- as.numeric(
  as.character(coords_root_600$GPSlongitude)
)

coords_root_600$GPSlatitude <- as.numeric(
  as.character(coords_root_600$GPSlatitude)
)


# Check whether sample IDs match
setequal(
  rownames(coords_root_600),
  rownames(biom_dm_root_600)
)


# Match coordinate order with UniFrac matrix
coords_root_600 <- coords_root_600[
  rownames(biom_dm_root_600),
  ,
  drop = FALSE
]

# Match UniFrac column order with row order
biom_dm_root_600 <- biom_dm_root_600[
  rownames(biom_dm_root_600),
  rownames(biom_dm_root_600)
]


# Verify sample order
identical(
  rownames(coords_root_600),
  rownames(biom_dm_root_600)
)


# Check GPS coordinates
summary(coords_root_600)

colSums(is.na(coords_root_600))


# Convert coordinates to matrix
coords_root_600 <- as.matrix(coords_root_600)


#-------------------------------------------------
# Compute Physical Distances (Haversine)
#-------------------------------------------------

# Calculate geographic distance in meters
phys_dist_matrix_root_600 <- distm(
  coords_root_600,
  fun = distHaversine
)

# Convert meters to kilometers
phys_dist_matrix_root_600 <- phys_dist_matrix_root_600 / 1000

phys_dist_root_600 <- as.dist(phys_dist_matrix_root_600)


#-------------------------------------------------
# Convert to Vectors for Correlation
#-------------------------------------------------

unifrac_vec_root_600 <- as.vector(
  as.dist(biom_dm_root_600)
)

phys_vec_root_600 <- as.vector(
  phys_dist_root_600
)


#-------------------------------------------------
# Check Geographic Distances (km)
#-------------------------------------------------

summary(phys_vec_root_600)

max(phys_vec_root_600)


#-------------------------------------------------
# Log-Transform Physical Distances
#-------------------------------------------------

phys_vec_log_root_600 <- log1p(phys_vec_root_600)


#-------------------------------------------------
# Spearman Correlation
#-------------------------------------------------

cor_result_root_600 <- cor.test(
  unifrac_vec_root_600,
  phys_vec_log_root_600,
  method = "spearman",
  exact = FALSE
)

print(cor_result_root_600)


#-------------------------------------------------
# Create Data Frame for Plotting
#-------------------------------------------------

df_dist_root_600 <- data.frame(
  UniFrac_Distance_root_600 = unifrac_vec_root_600,
  Physical_Distance_km_root_600 = phys_vec_root_600,
  Physical_Distance_log_root_600 = phys_vec_log_root_600
)

head(df_dist_root_600)


#------------------------------------------------------------------------------------------------------------
# Perform the Mantel test, Since you are analyzing geographic distance against weighted UniFrac distance
#---------------------------------------------------------------------------------------------------------
mantel_root_600 <- mantel(
  biom_dm_root_600,
  as.matrix(phys_dist_root_600),
  method = "spearman",
  permutations = 999
)

print(mantel_root_600)



#==========================================================
# root samples - 1400 GDD
#==========================================================

ps_root_1400 <- subset_samples(
  root_bac,
  Growing_Degree_Days == "1400"
)

root_tree_1400 <- phy_tree(ps_root_1400)

# Check if tree is binary
sprintf("Is tree binary: %s", is.binary(root_tree_1400))

# Resolve polytomies
phy_tree(ps_root_1400) <- multi2di(root_tree_1400)

sprintf("Is tree binary: %s", is.binary(phy_tree(ps_root_1400)))


# Convert phyloseq object to rbiom object
biom_root_1400 <- as_rbiom(ps_root_1400)
biom_root_1400

# View the biom file
glimpse(biom_root_1400, width = NULL)


# Calculate UniFrac distance between samples
biom_dm_root_1400 <- bdiv_distmat(
  biom_root_1400,
  bdiv = "UniFrac",
  weighted = TRUE,
  normalized = TRUE
)
biom_dm_root_1400


# Convert UniFrac distance to matrix
biom_dm_root_1400 <- as.matrix(biom_dm_root_1400)


#-------------------------------------------------
# Check sample IDs
#-------------------------------------------------

# Check if samples are the same
setequal(
  sample_names(ps_root_1400),
  rownames(biom_dm_root_1400)
)

# Check if sample order is the same
identical(
  sample_names(ps_root_1400),
  rownames(biom_dm_root_1400)
)

# Check UniFrac row and column names
identical(
  rownames(biom_dm_root_1400),
  colnames(biom_dm_root_1400)
)


#-------------------------------------------------
# Extract GPS coordinates
#-------------------------------------------------

# Longitude first, latitude second
coords_root_1400 <- sample_data(ps_root_1400)[,
                                                c("GPSlongitude", "GPSlatitude")
]

coords_root_1400 <- as.data.frame(coords_root_1400)

# Convert coordinates to numeric
coords_root_1400$GPSlongitude <- as.numeric(
  as.character(coords_root_1400$GPSlongitude)
)

coords_root_1400$GPSlatitude <- as.numeric(
  as.character(coords_root_1400$GPSlatitude)
)


# Check whether sample IDs match
setequal(
  rownames(coords_root_1400),
  rownames(biom_dm_root_1400)
)


# Match coordinate order with UniFrac matrix
coords_root_1400 <- coords_root_1400[
  rownames(biom_dm_root_1400),
  ,
  drop = FALSE
]

# Match UniFrac column order with row order
biom_dm_root_1400 <- biom_dm_root_1400[
  rownames(biom_dm_root_1400),
  rownames(biom_dm_root_1400)
]


# Verify sample order
identical(
  rownames(coords_root_1400),
  rownames(biom_dm_root_1400)
)


# Check GPS coordinates
summary(coords_root_1400)

colSums(is.na(coords_root_1400))


# Convert coordinates to matrix
coords_root_1400 <- as.matrix(coords_root_1400)


#-------------------------------------------------
# Compute Physical Distances (Haversine)
#-------------------------------------------------

# Calculate geographic distance in meters
phys_dist_matrix_root_1400 <- distm(
  coords_root_1400,
  fun = distHaversine
)

# Convert meters to kilometers
phys_dist_matrix_root_1400 <- phys_dist_matrix_root_1400 / 1000

phys_dist_root_1400 <- as.dist(phys_dist_matrix_root_1400)


#-------------------------------------------------
# Convert to Vectors for Correlation
#-------------------------------------------------

unifrac_vec_root_1400 <- as.vector(
  as.dist(biom_dm_root_1400)
)

phys_vec_root_1400 <- as.vector(
  phys_dist_root_1400
)


#-------------------------------------------------
# Check Geographic Distances (km)
#-------------------------------------------------

summary(phys_vec_root_1400)

max(phys_vec_root_1400)


#-------------------------------------------------
# Log-Transform Physical Distances
#-------------------------------------------------

phys_vec_log_root_1400 <- log1p(phys_vec_root_1400)


#-------------------------------------------------
# Spearman Correlation
#-------------------------------------------------

cor_result_root_1400 <- cor.test(
  unifrac_vec_root_1400,
  phys_vec_log_root_1400,
  method = "spearman",
  exact = FALSE
)

print(cor_result_root_1400)


#-------------------------------------------------
# Create Data Frame for Plotting
#-------------------------------------------------

df_dist_root_1400 <- data.frame(
  UniFrac_Distance_root_1400 = unifrac_vec_root_1400,
  Physical_Distance_km_root_1400 = phys_vec_root_1400,
  Physical_Distance_log_root_1400 = phys_vec_log_root_1400
)

head(df_dist_root_1400)


#------------------------------------------------------------------------------------------------------------
# Perform the Mantel test, Since you are analyzing geographic distance against weighted UniFrac distance
#---------------------------------------------------------------------------------------------------------
mantel_root_1400 <- mantel(
  biom_dm_root_1400,
  as.matrix(phys_dist_root_1400),
  method = "spearman",
  permutations = 999
)

print(mantel_root_1400)



#For Stalk samples


#==========================================================
# stalk samples - 600 GDD
#==========================================================

ps_stalk_600 <- subset_samples(
  stalk_bac,
  Growing_Degree_Days == "600"
)

stalk_tree_600 <- phy_tree(ps_stalk_600)

# Check if tree is binary
sprintf("Is tree binary: %s", is.binary(stalk_tree_600))

# Resolve polytomies
phy_tree(ps_stalk_600) <- multi2di(stalk_tree_600)

sprintf("Is tree binary: %s", is.binary(phy_tree(ps_stalk_600)))


# Convert phyloseq object to rbiom object
biom_stalk_600 <- as_rbiom(ps_stalk_600)
biom_stalk_600

# View the biom file
glimpse(biom_stalk_600, width = NULL)


# Calculate UniFrac distance between samples
biom_dm_stalk_600 <- bdiv_distmat(
  biom_stalk_600,
  bdiv = "UniFrac",
  weighted = TRUE,
  normalized = TRUE
)
biom_dm_stalk_600


# Convert UniFrac distance to matrix
biom_dm_stalk_600 <- as.matrix(biom_dm_stalk_600)


#-------------------------------------------------
# Check sample IDs
#-------------------------------------------------

# Check if samples are the same
setequal(
  sample_names(ps_stalk_600),
  rownames(biom_dm_stalk_600)
)

# Check if sample order is the same
identical(
  sample_names(ps_stalk_600),
  rownames(biom_dm_stalk_600)
)

# Check UniFrac row and column names
identical(
  rownames(biom_dm_stalk_600),
  colnames(biom_dm_stalk_600)
)


#-------------------------------------------------
# Extract GPS coordinates
#-------------------------------------------------

# Longitude first, latitude second
coords_stalk_600 <- sample_data(ps_stalk_600)[,
                                            c("GPSlongitude", "GPSlatitude")
]

coords_stalk_600 <- as.data.frame(coords_stalk_600)

# Convert coordinates to numeric
coords_stalk_600$GPSlongitude <- as.numeric(
  as.character(coords_stalk_600$GPSlongitude)
)

coords_stalk_600$GPSlatitude <- as.numeric(
  as.character(coords_stalk_600$GPSlatitude)
)


# Check whether sample IDs match
setequal(
  rownames(coords_stalk_600),
  rownames(biom_dm_stalk_600)
)


# Match coordinate order with UniFrac matrix
coords_stalk_600 <- coords_stalk_600[
  rownames(biom_dm_stalk_600),
  ,
  drop = FALSE
]

# Match UniFrac column order with row order
biom_dm_stalk_600 <- biom_dm_stalk_600[
  rownames(biom_dm_stalk_600),
  rownames(biom_dm_stalk_600)
]


# Verify sample order
identical(
  rownames(coords_stalk_600),
  rownames(biom_dm_stalk_600)
)


# Check GPS coordinates
summary(coords_stalk_600)

colSums(is.na(coords_stalk_600))


# Convert coordinates to matrix
coords_stalk_600 <- as.matrix(coords_stalk_600)


#-------------------------------------------------
# Compute Physical Distances (Haversine)
#-------------------------------------------------

# Calculate geographic distance in meters
phys_dist_matrix_stalk_600 <- distm(
  coords_stalk_600,
  fun = distHaversine
)

# Convert meters to kilometers
phys_dist_matrix_stalk_600 <- phys_dist_matrix_stalk_600 / 1000

phys_dist_stalk_600 <- as.dist(phys_dist_matrix_stalk_600)


#-------------------------------------------------
# Convert to Vectors for Correlation
#-------------------------------------------------

unifrac_vec_stalk_600 <- as.vector(
  as.dist(biom_dm_stalk_600)
)

phys_vec_stalk_600 <- as.vector(
  phys_dist_stalk_600
)


#-------------------------------------------------
# Check Geographic Distances (km)
#-------------------------------------------------

summary(phys_vec_stalk_600)

max(phys_vec_stalk_600)


#-------------------------------------------------
# Log-Transform Physical Distances
#-------------------------------------------------

phys_vec_log_stalk_600 <- log1p(phys_vec_stalk_600)


#-------------------------------------------------
# Spearman Correlation
#-------------------------------------------------

cor_result_stalk_600 <- cor.test(
  unifrac_vec_stalk_600,
  phys_vec_log_stalk_600,
  method = "spearman",
  exact = FALSE
)

print(cor_result_stalk_600)


#-------------------------------------------------
# Create Data Frame for Plotting
#-------------------------------------------------

df_dist_stalk_600 <- data.frame(
  UniFrac_Distance_stalk_600 = unifrac_vec_stalk_600,
  Physical_Distance_km_stalk_600 = phys_vec_stalk_600,
  Physical_Distance_log_stalk_600 = phys_vec_log_stalk_600
)

head(df_dist_stalk_600)


#------------------------------------------------------------------------------------------------------------
# Perform the Mantel test, Since you are analyzing geographic distance against weighted UniFrac distance
#---------------------------------------------------------------------------------------------------------
mantel_stalk_600 <- mantel(
  biom_dm_stalk_600,
  as.matrix(phys_dist_stalk_600),
  method = "spearman",
  permutations = 999
)

print(mantel_stalk_600)



#==========================================================
# stalk samples - 1400 GDD
#==========================================================

ps_stalk_1400 <- subset_samples(
  stalk_bac,
  Growing_Degree_Days == "1400"
)

stalk_tree_1400 <- phy_tree(ps_stalk_1400)

# Check if tree is binary
sprintf("Is tree binary: %s", is.binary(stalk_tree_1400))

# Resolve polytomies
phy_tree(ps_stalk_1400) <- multi2di(stalk_tree_1400)

sprintf("Is tree binary: %s", is.binary(phy_tree(ps_stalk_1400)))


# Convert phyloseq object to rbiom object
biom_stalk_1400 <- as_rbiom(ps_stalk_1400)
biom_stalk_1400

# View the biom file
glimpse(biom_stalk_1400, width = NULL)


# Calculate UniFrac distance between samples
biom_dm_stalk_1400 <- bdiv_distmat(
  biom_stalk_1400,
  bdiv = "UniFrac",
  weighted = TRUE,
  normalized = TRUE
)
biom_dm_stalk_1400


# Convert UniFrac distance to matrix
biom_dm_stalk_1400 <- as.matrix(biom_dm_stalk_1400)


#-------------------------------------------------
# Check sample IDs
#-------------------------------------------------

# Check if samples are the same
setequal(
  sample_names(ps_stalk_1400),
  rownames(biom_dm_stalk_1400)
)

# Check if sample order is the same
identical(
  sample_names(ps_stalk_1400),
  rownames(biom_dm_stalk_1400)
)

# Check UniFrac row and column names
identical(
  rownames(biom_dm_stalk_1400),
  colnames(biom_dm_stalk_1400)
)


#-------------------------------------------------
# Extract GPS coordinates
#-------------------------------------------------

# Longitude first, latitude second
coords_stalk_1400 <- sample_data(ps_stalk_1400)[,
                                              c("GPSlongitude", "GPSlatitude")
]

coords_stalk_1400 <- as.data.frame(coords_stalk_1400)

# Convert coordinates to numeric
coords_stalk_1400$GPSlongitude <- as.numeric(
  as.character(coords_stalk_1400$GPSlongitude)
)

coords_stalk_1400$GPSlatitude <- as.numeric(
  as.character(coords_stalk_1400$GPSlatitude)
)


# Check whether sample IDs match
setequal(
  rownames(coords_stalk_1400),
  rownames(biom_dm_stalk_1400)
)


# Match coordinate order with UniFrac matrix
coords_stalk_1400 <- coords_stalk_1400[
  rownames(biom_dm_stalk_1400),
  ,
  drop = FALSE
]

# Match UniFrac column order with row order
biom_dm_stalk_1400 <- biom_dm_stalk_1400[
  rownames(biom_dm_stalk_1400),
  rownames(biom_dm_stalk_1400)
]


# Verify sample order
identical(
  rownames(coords_stalk_1400),
  rownames(biom_dm_stalk_1400)
)


# Check GPS coordinates
summary(coords_stalk_1400)

colSums(is.na(coords_stalk_1400))


# Convert coordinates to matrix
coords_stalk_1400 <- as.matrix(coords_stalk_1400)


#-------------------------------------------------
# Compute Physical Distances (Haversine)
#-------------------------------------------------

# Calculate geographic distance in meters
phys_dist_matrix_stalk_1400 <- distm(
  coords_stalk_1400,
  fun = distHaversine
)

# Convert meters to kilometers
phys_dist_matrix_stalk_1400 <- phys_dist_matrix_stalk_1400 / 1000

phys_dist_stalk_1400 <- as.dist(phys_dist_matrix_stalk_1400)


#-------------------------------------------------
# Convert to Vectors for Correlation
#-------------------------------------------------

unifrac_vec_stalk_1400 <- as.vector(
  as.dist(biom_dm_stalk_1400)
)

phys_vec_stalk_1400 <- as.vector(
  phys_dist_stalk_1400
)


#-------------------------------------------------
# Check Geographic Distances (km)
#-------------------------------------------------

summary(phys_vec_stalk_1400)

max(phys_vec_stalk_1400)


#-------------------------------------------------
# Log-Transform Physical Distances
#-------------------------------------------------

phys_vec_log_stalk_1400 <- log1p(phys_vec_stalk_1400)


#-------------------------------------------------
# Spearman Correlation
#-------------------------------------------------

cor_result_stalk_1400 <- cor.test(
  unifrac_vec_stalk_1400,
  phys_vec_log_stalk_1400,
  method = "spearman",
  exact = FALSE
)

print(cor_result_stalk_1400)


#-------------------------------------------------
# Create Data Frame for Plotting
#-------------------------------------------------

df_dist_stalk_1400 <- data.frame(
  UniFrac_Distance_stalk_1400 = unifrac_vec_stalk_1400,
  Physical_Distance_km_stalk_1400 = phys_vec_stalk_1400,
  Physical_Distance_log_stalk_1400 = phys_vec_log_stalk_1400
)

head(df_dist_stalk_1400)


#------------------------------------------------------------------------------------------------------------
# Perform the Mantel test, Since you are analyzing geographic distance against weighted UniFrac distance
#---------------------------------------------------------------------------------------------------------
mantel_stalk_1400 <- mantel(
  biom_dm_stalk_1400,
  as.matrix(phys_dist_stalk_1400),
  method = "spearman",
  permutations = 999
)

print(mantel_stalk_1400)



#=========================================================
# Function to create sample-level distance-decay plots
#=========================================================

plot_distance_decay <- function(data,
                                x_var,
                                y_var,
                                title,
                                save_name,
                                mantel_result,
                                image_path = "C:/Aduragbemi/Manuscript/Review 8/IMAGES8",
                                pdf_path = "C:/Aduragbemi/Manuscript/Review 8/PDF8") {
  
  p <- ggplot(
    data,
    aes(
      x = .data[[x_var]],
      y = .data[[y_var]]
    )
  ) +
    
    theme_gray(base_size = 14) +
    
    geom_point(
      alpha = 0.4,
      size = 6
    ) +
    
    geom_smooth(
      method = "lm",
      se = FALSE,
      colour = "blue",
      linetype = "dashed"
    ) +
    
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = paste0(
        "Mantel's r = ",
        round(mantel_result$statistic, 3),
        ", p = ",
        round(mantel_result$signif, 3)
      ),
      hjust = -0.1,
      vjust = 1.5,
      size = 6
    ) +
    
    labs(
      title = title,
      x = "Geographic Distance (km)",
      y = "Weighted UniFrac Distance"
    ) +
    
    theme(
      
      plot.title = element_text(
        family = "serif",
        face = "bold",
        size = 15
      ),
      
      axis.title = element_text(
        family = "serif",
        face = "bold",
        size = 15,
        colour = "black"
      ),
      
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        size = 15,
        colour = "black"
      ),
      
      axis.text.y = element_text(
        angle = 0,
        hjust = 1,
        size = 15,
        colour = "black"
      ),
      
      legend.position = "bottom",
      
      legend.title = element_text(
        face = "bold",
        size = 15
      ),
      
      legend.text = element_text(
        size = 15,
        colour = "black"
      )
    )
  
  print(p)
  
  # Save PNG
  ggsave(
    paste0(save_name, ".png"),
    plot = p,
    path = image_path,
    dpi = 700,
    width = 10,
    height = 5,
    units = "in"
  )
  
  # Save PDF
  ggsave(
    paste0(save_name, ".pdf"),
    plot = p,
    path = pdf_path,
    width = 10,
    height = 5,
    units = "in"
  )
  
  return(p)
}


#=========================================================
# Soil samples - 600 GDD
#=========================================================

plot_soil_600 <- plot_distance_decay(
  data = df_dist_soil_600,
  x_var = "Physical_Distance_km_soil_600",
  y_var = "UniFrac_Distance_soil_600",
  title = "Soil microbiome - 600 GDD",
  save_name = "plot_soil_600",
  mantel_result = mantel_soil_600
)

plot_soil_600


#=========================================================
# Soil samples - 1400 GDD
#=========================================================

plot_soil_1400 <- plot_distance_decay(
  data = df_dist_soil_1400,
  x_var = "Physical_Distance_km_soil_1400",
  y_var = "UniFrac_Distance_soil_1400",
  title = "Soil microbiome - 1400 GDD",
  save_name = "plot_soil_1400",
  mantel_result = mantel_soil_1400
)

plot_soil_1400


#=========================================================
# Rhizosphere samples - 600 GDD
#=========================================================

plot_rhizo_600 <- plot_distance_decay(
  data = df_dist_rhizo_600,
  x_var = "Physical_Distance_km_rhizo_600",
  y_var = "UniFrac_Distance_rhizo_600",
  title = "Rhizosphere microbiome - 600 GDD",
  save_name = "plot_rhizo_600",
  mantel_result = mantel_rhizo_600
)

plot_rhizo_600


#=========================================================
# Rhizosphere samples - 1400 GDD
#=========================================================

plot_rhizo_1400 <- plot_distance_decay(
  data = df_dist_rhizo_1400,
  x_var = "Physical_Distance_km_rhizo_1400",
  y_var = "UniFrac_Distance_rhizo_1400",
  title = "Rhizosphere microbiome - 1400 GDD",
  save_name = "plot_rhizo_1400",
  mantel_result = mantel_rhizo_1400
)

plot_rhizo_1400


#=========================================================
# Root samples - 600 GDD
#=========================================================

plot_root_600 <- plot_distance_decay(
  data = df_dist_root_600,
  x_var = "Physical_Distance_km_root_600",
  y_var = "UniFrac_Distance_root_600",
  title = "Root microbiome - 600 GDD",
  save_name = "plot_root_600",
  mantel_result = mantel_root_600
)

plot_root_600


#=========================================================
# Root samples - 1400 GDD
#=========================================================

plot_root_1400 <- plot_distance_decay(
  data = df_dist_root_1400,
  x_var = "Physical_Distance_km_root_1400",
  y_var = "UniFrac_Distance_root_1400",
  title = "Root microbiome - 1400 GDD",
  save_name = "plot_root_1400",
  mantel_result = mantel_root_1400
)

plot_root_1400


#=========================================================
# Stalk samples - 600 GDD
#=========================================================

plot_stalk_600 <- plot_distance_decay(
  data = df_dist_stalk_600,
  x_var = "Physical_Distance_km_stalk_600",
  y_var = "UniFrac_Distance_stalk_600",
  title = "Stalk microbiome - 600 GDD",
  save_name = "plot_stalk_600",
  mantel_result = mantel_stalk_600
)

plot_stalk_600


#=========================================================
# Stalk samples - 1400 GDD
#=========================================================

plot_stalk_1400 <- plot_distance_decay(
  data = df_dist_stalk_1400,
  x_var = "Physical_Distance_km_stalk_1400",
  y_var = "UniFrac_Distance_stalk_1400",
  title = "Stalk microbiome - 1400 GDD",
  save_name = "plot_stalk_1400",
  mantel_result = mantel_stalk_1400
)

plot_stalk_1400
