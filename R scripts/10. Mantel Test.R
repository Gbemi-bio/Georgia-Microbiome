# Load necessary libraries
library(phyloseq)
library(ape)
library(vegan)
library(geosphere)
library(rbiom)

set.seed(1234)


#For Soil Sample
#read in the data
soil <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/soil_bac_sub_no_bad_Filtered_final_non_normalized.rds")


# Convert phyloseq object to rbiom object
soil_physeq_biom <- as_rbiom(soil)
soil_physeq_biom

#view the biom file
glimpse(soil_physeq_biom, width = NULL)

# Calculate weighted UniFrac distance
dist_soil <- bdiv_distmat(soil_physeq_biom, 'unifrac', weighted = TRUE)
dist_soil

# Print the distance matrix
print(dist_soil)

# Extract the sample data
sample_data_soil <- sample_data(soil)

# Convert sample_data to a data frame (if it is not already)
sample_data_soil_df <- as.data.frame(sample_data_soil)

#for geography
# Subset the desired columns for geography
subset_sample_data_soil <- sample_data_soil_df[, c("GPSlatitude", "GPSlongitude")]

# Print the subset of the sample data
print(subset_sample_data_soil)

#geographic data frame - haversine distance 
subset_sample_data_soil = distm(subset_sample_data_soil, fun = distHaversine)
dist.subset_sample_data_soil = as.dist(subset_sample_data_soil)

#for environmental factors
# Subset the desired columns for environmental variables
soil_subset_sample_data_env <- sample_data_soil_df[, c("Growing_Degree_Days", "precipIntensity", "precipProbability", "temperatureHigh",
                                                       "temperatureLow", "dewPoint", "humidity", "windSpeed", "uvIndex")]

# Print the subset of  environmental factors
print(soil_subset_sample_data_env)


# Inspect the data types of each column
str(soil_subset_sample_data_env)

# Convert non-numeric columns to numeric (if necessary)
soil_subset_sample_data_env <- as.data.frame(lapply(soil_subset_sample_data_env, function(x) as.numeric(as.character(x))))

# Verify the conversion
str(soil_subset_sample_data_env)

#scale data environmental data frame 
scale.env_soil = scale(soil_subset_sample_data_env, center = TRUE, scale = TRUE)

#create distance matrix of scaled data
dist.env_soil = dist(scale.env_soil, method = "euclidean")


#for soil properties
# Subset the desired columns for soil properties
subset_sample_data_soil_prop <- sample_data_soil_df[, c("Lime_Buffer_Capacity", "pH", "Water_pH", "Base_Saturation_perc",
                                                        "Cation_Exchange_Capacity", "Ca_ppm", "Cd_ppm", "Cr_ppm", "Cu_ppm",
                                                        "Fe_ppm", "K_ppm", "Mg_ppm", "Mn_ppm", "Mo_ppm", "Na_ppm", "Ni_ppm",
                                                        "P_ppm", "Pb_ppm", "Zn_ppm", "Organic_Matter_perc", "N_perc")]

# Inspect the data types of each column
str(subset_sample_data_soil_prop)

# Convert non-numeric columns to numeric
subset_sample_data_soil_prop <- as.data.frame(lapply(subset_sample_data_soil_prop, function(x) as.numeric(as.character(x))))

# Verify the conversion
str(subset_sample_data_soil_prop)

#scale data soil properties data frame 
scale.prop_soil = scale(subset_sample_data_soil_prop, center = TRUE, scale = TRUE)

#create distance matrix of scaled data
dist.prop_soil = dist(scale.prop_soil, method = "euclidean")


#now run mantel test for species abundance vs geography
abund_geo_soil  = mantel(dist_soil, dist.subset_sample_data_soil, method = "spearman", permutations = 999, na.rm = TRUE)
abund_geo_soil

#run mantel test for species abundance vs soil chemistry 
abund_prop_soil = mantel(dist_soil, dist.prop_soil, method = "spearman", permutations = 999, na.rm = TRUE)
abund_prop_soil

#run mantel test for species abundance vs environmental factors
abund_env_soil = mantel(dist_soil, dist.env_soil, method = "spearman", permutations = 999, na.rm = TRUE)
abund_env_soil







#For Rhizosphere Sample
#read in the data
rhizo <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/rhizo_bac_sub_no_bad_Filtered_final_non_normalized.rds")


# Convert phyloseq object to rbiom object
rhizo_physeq_biom <- as_rbiom(rhizo)
rhizo_physeq_biom

#view the biom file
glimpse(rhizo_physeq_biom, width = NULL)

# Calculate weighted UniFrac distance
dist_rhizo <- bdiv_distmat(rhizo_physeq_biom, 'unifrac', weighted = TRUE)
dist_rhizo

# Print the distance matrix
print(dist_rhizo)

# Extract the sample data
sample_data_rhizo <- sample_data(rhizo)

# Convert sample_data to a data frame (if it is not already)
sample_data_rhizo_df <- as.data.frame(sample_data_rhizo)

#for geography
# Subset the desired columns for geography
subset_sample_data_rhizo <- sample_data_rhizo_df[, c("GPSlatitude", "GPSlongitude")]

# Print the subset of the sample data
print(subset_sample_data_rhizo)

#geographic data frame - haversine distance 
subset_sample_data_rhizo = distm(subset_sample_data_rhizo, fun = distHaversine)
dist.subset_sample_data_rhizo = as.dist(subset_sample_data_rhizo)

#for environmental factors
# Subset the desired columns for environmental variables
rhizo_subset_sample_data_env <- sample_data_rhizo_df[, c("Growing_Degree_Days", "precipIntensity", "precipProbability", "temperatureHigh",
                                                       "temperatureLow", "dewPoint", "humidity", "windSpeed", "uvIndex")]

# Print the subset of  environmental factors
print(rhizo_subset_sample_data_env)


# Inspect the data types of each column
str(rhizo_subset_sample_data_env)

# Convert non-numeric columns to numeric (if necessary)
rhizo_subset_sample_data_env <- as.data.frame(lapply(rhizo_subset_sample_data_env, function(x) as.numeric(as.character(x))))

# Verify the conversion
str(rhizo_subset_sample_data_env)

#scale data environmental data frame 
scale.env_rhizo = scale(rhizo_subset_sample_data_env, center = TRUE, scale = TRUE)

#create distance matrix of scaled data
dist.env_rhizo = dist(scale.env_rhizo, method = "euclidean")


#for rhizo properties
# Subset the desired columns for rhizo properties
subset_sample_data_rhizo_prop <- sample_data_rhizo_df[, c("Lime_Buffer_Capacity", "pH", "Water_pH", "Base_Saturation_perc",
                                                        "Cation_Exchange_Capacity", "Ca_ppm", "Cd_ppm", "Cr_ppm", "Cu_ppm",
                                                        "Fe_ppm", "K_ppm", "Mg_ppm", "Mn_ppm", "Mo_ppm", "Na_ppm", "Ni_ppm",
                                                        "P_ppm", "Pb_ppm", "Zn_ppm", "Organic_Matter_perc", "N_perc")]

# Inspect the data types of each column
str(subset_sample_data_rhizo_prop)

# Convert non-numeric columns to numeric
subset_sample_data_rhizo_prop <- as.data.frame(lapply(subset_sample_data_rhizo_prop, function(x) as.numeric(as.character(x))))

# Verify the conversion
str(subset_sample_data_rhizo_prop)

#scale data rhizo properties data frame 
scale.prop_rhizo = scale(subset_sample_data_rhizo_prop, center = TRUE, scale = TRUE)

#create distance matrix of scaled data
dist.prop_rhizo = dist(scale.prop_rhizo, method = "euclidean")


#now run mantel test for species abundance vs geography
abund_geo_rhizo  = mantel(dist_rhizo, dist.subset_sample_data_rhizo, method = "spearman", permutations = 999, na.rm = TRUE)
abund_geo_rhizo

#run mantel test for species abundance vs soil chemistry 
abund_prop_rhizo = mantel(dist_rhizo, dist.prop_rhizo, method = "spearman", permutations = 999, na.rm = TRUE)
abund_prop_rhizo

#run mantel test for species abundance vs environmental factors
abund_env_rhizo = mantel(dist_rhizo, dist.env_rhizo, method = "spearman", permutations = 999, na.rm = TRUE)
abund_env_rhizo



#For Root Sample
#read in the data
root <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/root_bac_sub_no_bad_Filtered_final_non_normalized.rds")


# Convert phyloseq object to rbiom object
root_physeq_biom <- as_rbiom(root)
root_physeq_biom

#view the biom file
glimpse(root_physeq_biom, width = NULL)

# Calculate weighted UniFrac distance
dist_root <- bdiv_distmat(root_physeq_biom, 'unifrac', weighted = TRUE)
dist_root

# Print the distance matrix
print(dist_root)

# Extract the sample data
sample_data_root <- sample_data(root)

# Convert sample_data to a data frame (if it is not already)
sample_data_root_df <- as.data.frame(sample_data_root)

#for geography
# Subset the desired columns for geography
subset_sample_data_root <- sample_data_root_df[, c("GPSlatitude", "GPSlongitude")]

# Print the subset of the sample data
print(subset_sample_data_root)

#geographic data frame - haversine distance 
subset_sample_data_root = distm(subset_sample_data_root, fun = distHaversine)
dist.subset_sample_data_root = as.dist(subset_sample_data_root)

#for environmental factors
# Subset the desired columns for environmental variables
root_subset_sample_data_env <- sample_data_root_df[, c("Growing_Degree_Days", "precipIntensity", "precipProbability", "temperatureHigh",
                                                         "temperatureLow", "dewPoint", "humidity", "windSpeed", "uvIndex")]

# Print the subset of  environmental factors
print(root_subset_sample_data_env)


# Inspect the data types of each column
str(root_subset_sample_data_env)

# Convert non-numeric columns to numeric (if necessary)
root_subset_sample_data_env <- as.data.frame(lapply(root_subset_sample_data_env, function(x) as.numeric(as.character(x))))

# Verify the conversion
str(root_subset_sample_data_env)

#scale data environmental data frame 
scale.env_root = scale(root_subset_sample_data_env, center = TRUE, scale = TRUE)

#create distance matrix of scaled data
dist.env_root = dist(scale.env_root, method = "euclidean")


#for root properties
# Subset the desired columns for root properties
subset_sample_data_root_prop <- sample_data_root_df[, c("Lime_Buffer_Capacity", "pH", "Water_pH", "Base_Saturation_perc",
                                                          "Cation_Exchange_Capacity", "Ca_ppm", "Cd_ppm", "Cr_ppm", "Cu_ppm",
                                                          "Fe_ppm", "K_ppm", "Mg_ppm", "Mn_ppm", "Mo_ppm", "Na_ppm", "Ni_ppm",
                                                          "P_ppm", "Pb_ppm", "Zn_ppm", "Organic_Matter_perc", "N_perc")]

# Inspect the data types of each column
str(subset_sample_data_root_prop)

# Convert non-numeric columns to numeric
subset_sample_data_root_prop <- as.data.frame(lapply(subset_sample_data_root_prop, function(x) as.numeric(as.character(x))))

# Verify the conversion
str(subset_sample_data_root_prop)

#scale data root properties data frame 
scale.prop_root = scale(subset_sample_data_root_prop, center = TRUE, scale = TRUE)

#create distance matrix of scaled data
dist.prop_root = dist(scale.prop_root, method = "euclidean")


#now run mantel test for species abundance vs geography
abund_geo_root  = mantel(dist_root, dist.subset_sample_data_root, method = "spearman", permutations = 999, na.rm = TRUE)
abund_geo_root

#run mantel test for species abundance vs soil chemistry 
abund_prop_root = mantel(dist_root, dist.prop_root, method = "spearman", permutations = 999, na.rm = TRUE)
abund_prop_root

#run mantel test for species abundance vs environmental factors
abund_env_root = mantel(dist_root, dist.env_root, method = "spearman", permutations = 999, na.rm = TRUE)
abund_env_root






#For Stalk Sample
#read in the data
stalk <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/stalk_bac_sub_no_bad_Filtered_final_non_normalized.rds")


# Convert phyloseq object to rbiom object
stalk_physeq_biom <- as_rbiom(stalk)
stalk_physeq_biom

#view the biom file
glimpse(stalk_physeq_biom, width = NULL)

# Calculate weighted UniFrac distance
dist_stalk <- bdiv_distmat(stalk_physeq_biom, 'unifrac', weighted = TRUE)
dist_stalk

# Print the distance matrix
print(dist_stalk)

# Extract the sample data
sample_data_stalk <- sample_data(stalk)

# Convert sample_data to a data frame (if it is not already)
sample_data_stalk_df <- as.data.frame(sample_data_stalk)

#for geography
# Subset the desired columns for geography
subset_sample_data_stalk <- sample_data_stalk_df[, c("GPSlatitude", "GPSlongitude")]

# Print the subset of the sample data
print(subset_sample_data_stalk)

#geographic data frame - haversine distance 
subset_sample_data_stalk = distm(subset_sample_data_stalk, fun = distHaversine)
dist.subset_sample_data_stalk = as.dist(subset_sample_data_stalk)

#for environmental factors
# Subset the desired columns for environmental variables
stalk_subset_sample_data_env <- sample_data_stalk_df[, c("Growing_Degree_Days", "precipIntensity", "precipProbability", "temperatureHigh",
                                                       "temperatureLow", "dewPoint", "humidity", "windSpeed", "uvIndex")]

# Print the subset of  environmental factors
print(stalk_subset_sample_data_env)


# Inspect the data types of each column
str(stalk_subset_sample_data_env)

# Convert non-numeric columns to numeric (if necessary)
stalk_subset_sample_data_env <- as.data.frame(lapply(stalk_subset_sample_data_env, function(x) as.numeric(as.character(x))))

# Verify the conversion
str(stalk_subset_sample_data_env)

#scale data environmental data frame 
scale.env_stalk = scale(stalk_subset_sample_data_env, center = TRUE, scale = TRUE)

#create distance matrix of scaled data
dist.env_stalk = dist(scale.env_stalk, method = "euclidean")


#for stalk properties
# Subset the desired columns for stalk properties
subset_sample_data_stalk_prop <- sample_data_stalk_df[, c("Lime_Buffer_Capacity", "pH", "Water_pH", "Base_Saturation_perc",
                                                        "Cation_Exchange_Capacity", "Ca_ppm", "Cd_ppm", "Cr_ppm", "Cu_ppm",
                                                        "Fe_ppm", "K_ppm", "Mg_ppm", "Mn_ppm", "Mo_ppm", "Na_ppm", "Ni_ppm",
                                                        "P_ppm", "Pb_ppm", "Zn_ppm", "Organic_Matter_perc", "N_perc")]

# Inspect the data types of each column
str(subset_sample_data_stalk_prop)

# Convert non-numeric columns to numeric
subset_sample_data_stalk_prop <- as.data.frame(lapply(subset_sample_data_stalk_prop, function(x) as.numeric(as.character(x))))

# Verify the conversion
str(subset_sample_data_stalk_prop)

#scale data stalk properties data frame 
scale.prop_stalk = scale(subset_sample_data_stalk_prop, center = TRUE, scale = TRUE)

#create distance matrix of scaled data
dist.prop_stalk = dist(scale.prop_stalk, method = "euclidean")


#now run mantel test for species abundance vs geography
abund_geo_stalk  = mantel(dist_stalk, dist.subset_sample_data_stalk, method = "spearman", permutations = 999, na.rm = TRUE)
abund_geo_stalk

#run mantel test for species abundance vs soil chemistry 
abund_prop_stalk = mantel(dist_stalk, dist.prop_stalk, method = "spearman", permutations = 999, na.rm = TRUE)
abund_prop_stalk

#run mantel test for species abundance vs environmental factors
abund_env_stalk = mantel(dist_stalk, dist.env_stalk, method = "spearman", permutations = 999, na.rm = TRUE)
abund_env_stalk




