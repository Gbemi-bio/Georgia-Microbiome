#Load required libraries
library(phyloseq)
library(ggplot2)
library(picante)

set.seed(1234)

# ================================
# read in the non-normalized data
# ================================
bac.css.non_norm_root <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/root_bac_sub_no_bad_Filtered_final_non_normalized.rds")
bac.css.non_norm_stalk <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/stalk_bac_sub_no_bad_Filtered_final_non_normalized.rds")
bac.css.non_norm_soil <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/soil_bac_sub_no_bad_Filtered_final_non_normalized.rds")
bac.css.non_norm_rhizo <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/rhizo_bac_sub_no_bad_Filtered_final_non_normalized.rds")

# ===============================================================================
# Rarefy to 90% of sample depth
# Rarefaction depth chosen is the 90% of the minimum sample depth in the dataset
# ===============================================================================

bac.css.non_norm_root_rare <- rarefy_even_depth(bac.css.non_norm_root, sample.size=1000, replace=F)
bac.css.non_norm_root_rare

bac.css.non_norm_soil_rare <- rarefy_even_depth(bac.css.non_norm_soil, sample.size=2000, replace=F)
bac.css.non_norm_soil_rare

bac.css.non_norm_stalk_rare <- rarefy_even_depth(bac.css.non_norm_stalk, sample.size=500, replace=F)
bac.css.non_norm_stalk_rare

bac.css.non_norm_rhizo_rare <- rarefy_even_depth(bac.css.non_norm_rhizo, sample.size=3000, replace=F)
bac.css.non_norm_rhizo_rare


bac.css.non_norm_soil_rare
bac.css.non_norm_stalk_rare
bac.css.non_norm_rhizo_rare
bac.css.non_norm_root_rare


# ======================================================================================
# Perform Alpha diversity of environmental factors,and geography on various sample types
# ======================================================================================

#For Soil sample

soil_continuous_alpha <- data.frame(bac.css.non_norm_soil_rare@sam_data)
soil_continuous_alpha
soil_continuous_alpha_1 <- estimate_richness(bac.css.non_norm_soil_rare, measures = c("Observed", "Shannon"))
soil_continuous_alpha_1

# Calculate Faith's Phylogenetic Diversity (PD)
pd_soil_otu <- as.data.frame(bac.css.non_norm_soil_rare@otu_table)

pd_soil_tree <- bac.css.non_norm_soil_rare@phy_tree

#check if the tree is rooted
bac.css.non_norm_soil_rare@phy_tree

# it is a rooted tree
faith_pd_soil <- pd(t(pd_soil_otu), pd_soil_tree,include.root=T)
faith_pd_soil
soil_continuous_alpha_1 <- cbind(soil_continuous_alpha, soil_continuous_alpha_1,faith_pd_soil)
soil_continuous_alpha_1 <- data.frame(soil_continuous_alpha_1)
soil_continuous_alpha_1


#for high temperature
glm_shannon_soil_hightemperature = glm(Shannon ~ temperatureHigh, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_hightemperature)

glm_shannon_soil_hightemperature = glm(PD ~ temperatureHigh, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_hightemperature)

#for low temperature
glm_shannon_soil_lowtemperature = glm(Shannon ~ temperatureLow, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_lowtemperature)

glm_shannon_soil_lowtemperature = glm(PD ~ temperatureLow, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_lowtemperature)


#for precipIntensity
glm_shannon_soil_precipIntensity = glm(Shannon ~ precipIntensity, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_precipIntensity)

glm_shannon_soil_precipIntensity = glm(PD ~ precipIntensity, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_precipIntensity)

#for windSpeed
glm_shannon_soil_windSpeed = glm(Shannon ~ windSpeed, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_windSpeed)

glm_shannon_soil_windSpeed = glm(PD ~ windSpeed, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_windSpeed)


#for humidity
glm_shannon_soil_humidity = glm(Shannon ~ humidity, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_humidity)

glm_shannon_soil_humidity = glm(PD ~ humidity, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_humidity)


#for uvindex
glm_shannon_soil_uvIndex = glm(Shannon ~ uvIndex, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_uvIndex)

glm_shannon_soil_uvIndex = glm(PD ~ uvIndex, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_uvIndex)

#for GPSlatitude
glm_shannon_soil_GPSlatitude = glm(Shannon ~ GPSlatitude, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_GPSlatitude)

glm_shannon_soil_GPSlatitude = glm(PD ~ GPSlatitude, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_GPSlatitude)

#for dewpoint
glm_shannon_soil_dewPoint = glm(Shannon ~ dewPoint, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_dewPoint)

glm_shannon_soil_dewPoint = glm(PD ~ dewPoint, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_dewPoint)


#for GPSlongitude
glm_shannon_soil_GPSlongitude = glm(Shannon ~ GPSlongitude, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_GPSlongitude)

glm_shannon_soil_GPSlongitude = glm(PD ~ GPSlongitude, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_GPSlongitude)


#For Rhizosphere sample

rhizo_continuous_alpha <- data.frame(bac.css.non_norm_rhizo_rare@sam_data)
rhizo_continuous_alpha
rhizo_continuous_alpha_1 <- estimate_richness(bac.css.non_norm_rhizo_rare, measures = c("Observed", "Shannon"))
rhizo_continuous_alpha_1

# Calculate Faith's Phylogenetic Diversity (PD)
pd_rhizo_otu <- as.data.frame(bac.css.non_norm_rhizo_rare@otu_table)

pd_rhizo_tree <- bac.css.non_norm_rhizo_rare@phy_tree

#check if the tree is rooted

bac.css.non_norm_rhizo_rare@phy_tree

# it is a rooted tree
faith_pd_rhizo <- pd(t(pd_rhizo_otu), pd_rhizo_tree,include.root=T)
faith_pd_rhizo
rhizo_continuous_alpha_1 <- cbind(rhizo_continuous_alpha, rhizo_continuous_alpha_1,faith_pd_rhizo)
rhizo_continuous_alpha_1 <- data.frame(rhizo_continuous_alpha_1)
rhizo_continuous_alpha_1


#for high temperature
glm_shannon_rhizo_hightemperature = glm(Shannon ~ temperatureHigh, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_hightemperature)

glm_shannon_rhizo_hightemperature = glm(PD ~ temperatureHigh, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_hightemperature)

#for low temperature
glm_shannon_rhizo_lowtemperature = glm(Shannon ~ temperatureLow, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_lowtemperature)

glm_shannon_rhizo_lowtemperature = glm(PD ~ temperatureLow, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_lowtemperature)


#for precipIntensity
glm_shannon_rhizo_precipIntensity = glm(Shannon ~ precipIntensity, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_precipIntensity)

glm_shannon_rhizo_precipIntensity = glm(PD ~ precipIntensity, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_precipIntensity)

#for windSpeed
glm_shannon_rhizo_windSpeed = glm(Shannon ~ windSpeed, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_windSpeed)

glm_shannon_rhizo_windSpeed = glm(PD ~ windSpeed, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_windSpeed)


#for humidity
glm_shannon_rhizo_humidity = glm(Shannon ~ humidity, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_humidity)

glm_shannon_rhizo_humidity = glm(PD ~ humidity, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_humidity)


#for uvindex
glm_shannon_rhizo_uvIndex = glm(Shannon ~ uvIndex, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_uvIndex)

glm_shannon_rhizo_uvIndex = glm(PD ~ uvIndex, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_uvIndex)

#for GPSlatitude
glm_shannon_rhizo_GPSlatitude = glm(Shannon ~ GPSlatitude, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_GPSlatitude)

glm_shannon_rhizo_GPSlatitude = glm(PD ~ GPSlatitude, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_GPSlatitude)

#for dewpoint
glm_shannon_rhizo_dewPoint = glm(Shannon ~ dewPoint, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_dewPoint)

glm_shannon_rhizo_dewPoint = glm(PD ~ dewPoint, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_dewPoint)


#for GPSlongitude
glm_shannon_rhizo_GPSlongitude = glm(Shannon ~ GPSlongitude, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_GPSlongitude)

glm_shannon_rhizo_GPSlongitude = glm(PD ~ GPSlongitude, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_GPSlongitude)


#For Root sample

root_continuous_alpha <- data.frame(bac.css.non_norm_root_rare@sam_data)
root_continuous_alpha
root_continuous_alpha_1 <- estimate_richness(bac.css.non_norm_root_rare, measures = c("Observed", "Shannon"))
root_continuous_alpha_1

# Calculate Faith's Phylogenetic Diversity (PD)
pd_root_otu <- as.data.frame(bac.css.non_norm_root_rare@otu_table)

pd_root_tree <- bac.css.non_norm_root_rare@phy_tree

#check if the tree is rooted

bac.css.non_norm_root_rare@phy_tree

# it is a rooted tree
faith_pd_root <- pd(t(pd_root_otu), pd_root_tree,include.root=T)
faith_pd_root
root_continuous_alpha_1 <- cbind(root_continuous_alpha, root_continuous_alpha_1,faith_pd_root)
root_continuous_alpha_1 <- data.frame(root_continuous_alpha_1)
root_continuous_alpha_1


#for high temperature
glm_shannon_root_hightemperature = glm(Shannon ~ temperatureHigh, data=root_continuous_alpha_1)
summary(glm_shannon_root_hightemperature)

glm_shannon_root_hightemperature = glm(PD ~ temperatureHigh, data=root_continuous_alpha_1)
summary(glm_shannon_root_hightemperature)

#for low temperature
glm_shannon_root_lowtemperature = glm(Shannon ~ temperatureLow, data=root_continuous_alpha_1)
summary(glm_shannon_root_lowtemperature)

glm_shannon_root_lowtemperature = glm(PD ~ temperatureLow, data=root_continuous_alpha_1)
summary(glm_shannon_root_lowtemperature)


#for precipIntensity
glm_shannon_root_precipIntensity = glm(Shannon ~ precipIntensity, data=root_continuous_alpha_1)
summary(glm_shannon_root_precipIntensity)

glm_shannon_root_precipIntensity = glm(PD ~ precipIntensity, data=root_continuous_alpha_1)
summary(glm_shannon_root_precipIntensity)

#for windSpeed
glm_shannon_root_windSpeed = glm(Shannon ~ windSpeed, data=root_continuous_alpha_1)
summary(glm_shannon_root_windSpeed)

glm_shannon_root_windSpeed = glm(PD ~ windSpeed, data=root_continuous_alpha_1)
summary(glm_shannon_root_windSpeed)


#for humidity
glm_shannon_root_humidity = glm(Shannon ~ humidity, data=root_continuous_alpha_1)
summary(glm_shannon_root_humidity)

glm_shannon_root_humidity = glm(PD ~ humidity, data=root_continuous_alpha_1)
summary(glm_shannon_root_humidity)


#for uvindex
glm_shannon_root_uvIndex = glm(Shannon ~ uvIndex, data=root_continuous_alpha_1)
summary(glm_shannon_root_uvIndex)

glm_shannon_root_uvIndex = glm(PD ~ uvIndex, data=root_continuous_alpha_1)
summary(glm_shannon_root_uvIndex)

#for GPSlatitude
glm_shannon_root_GPSlatitude = glm(Shannon ~ GPSlatitude, data=root_continuous_alpha_1)
summary(glm_shannon_root_GPSlatitude)

glm_shannon_root_GPSlatitude = glm(PD ~ GPSlatitude, data=root_continuous_alpha_1)
summary(glm_shannon_root_GPSlatitude)

#for dewpoint
glm_shannon_root_dewPoint = glm(Shannon ~ dewPoint, data=root_continuous_alpha_1)
summary(glm_shannon_root_dewPoint)

glm_shannon_root_dewPoint = glm(PD ~ dewPoint, data=root_continuous_alpha_1)
summary(glm_shannon_root_dewPoint)


#for GPSlongitude
glm_shannon_root_GPSlongitude = glm(Shannon ~ GPSlongitude, data=root_continuous_alpha_1)
summary(glm_shannon_root_GPSlongitude)

glm_shannon_root_GPSlongitude = glm(PD ~ GPSlongitude, data=root_continuous_alpha_1)
summary(glm_shannon_root_GPSlongitude)




#For Stalk sample

stalk_continuous_alpha <- data.frame(bac.css.non_norm_stalk_rare@sam_data)
stalk_continuous_alpha
stalk_continuous_alpha_1 <- estimate_richness(bac.css.non_norm_stalk_rare, measures = c("Observed", "Shannon"))
stalk_continuous_alpha_1

# Calculate Faith's Phylogenetic Diversity (PD)
pd_stalk_otu <- as.data.frame(bac.css.non_norm_stalk_rare@otu_table)

pd_stalk_tree <- bac.css.non_norm_stalk_rare@phy_tree

#check if the tree is stalked

bac.css.non_norm_stalk_rare@phy_tree

# it is a stalked tree
faith_pd_stalk <- pd(t(pd_stalk_otu), pd_stalk_tree,include.root=T)
faith_pd_stalk
stalk_continuous_alpha_1 <- cbind(stalk_continuous_alpha, stalk_continuous_alpha_1,faith_pd_stalk)
stalk_continuous_alpha_1 <- data.frame(stalk_continuous_alpha_1)
stalk_continuous_alpha_1


#for high temperature
glm_shannon_stalk_hightemperature = glm(Shannon ~ temperatureHigh, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_hightemperature)

glm_shannon_stalk_hightemperature = glm(PD ~ temperatureHigh, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_hightemperature)

#for low temperature
glm_shannon_stalk_lowtemperature = glm(Shannon ~ temperatureLow, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_lowtemperature)

glm_shannon_stalk_lowtemperature = glm(PD ~ temperatureLow, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_lowtemperature)


#for precipIntensity
glm_shannon_stalk_precipIntensity = glm(Shannon ~ precipIntensity, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_precipIntensity)

glm_shannon_stalk_precipIntensity = glm(PD ~ precipIntensity, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_precipIntensity)

#for windSpeed
glm_shannon_stalk_windSpeed = glm(Shannon ~ windSpeed, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_windSpeed)

glm_shannon_stalk_windSpeed = glm(PD ~ windSpeed, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_windSpeed)


#for humidity
glm_shannon_stalk_humidity = glm(Shannon ~ humidity, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_humidity)

glm_shannon_stalk_humidity = glm(PD ~ humidity, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_humidity)


#for uvindex
glm_shannon_stalk_uvIndex = glm(Shannon ~ uvIndex, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_uvIndex)

glm_shannon_stalk_uvIndex = glm(PD ~ uvIndex, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_uvIndex)

#for GPSlatitude
glm_shannon_stalk_GPSlatitude = glm(Shannon ~ GPSlatitude, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_GPSlatitude)

glm_shannon_stalk_GPSlatitude = glm(PD ~ GPSlatitude, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_GPSlatitude)

#for dewpoint
glm_shannon_stalk_dewPoint = glm(Shannon ~ dewPoint, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_dewPoint)

glm_shannon_stalk_dewPoint = glm(PD ~ dewPoint, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_dewPoint)


#for GPSlongitude
glm_shannon_stalk_GPSlongitude = glm(Shannon ~ GPSlongitude, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_GPSlongitude)

glm_shannon_stalk_GPSlongitude = glm(PD ~ GPSlongitude, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_GPSlongitude)


# =========================================================
# alpha diversity of Soil chemistry on various sample types
# =========================================================

#For soil sample

soil_continuous_alpha <- data.frame(bac.css.non_norm_soil_rare@sam_data)
soil_continuous_alpha
soil_continuous_alpha_1 <- estimate_richness(bac.css.non_norm_soil_rare, measures = c("Observed", "Shannon"))
soil_continuous_alpha_1

# Calculate Faith's Phylogenetic Diversity (PD)
pd_soil_otu <- as.data.frame(bac.css.non_norm_soil_rare@otu_table)

pd_soil_tree <- bac.css.non_norm_soil_rare@phy_tree

#check if the tree is rooted
bac.css.non_norm_soil_rare@phy_tree

# it is a rooted tree
faith_pd_soil <- pd(t(pd_soil_otu), pd_soil_tree,include.root=T)
faith_pd_soil
soil_continuous_alpha_1 <- cbind(soil_continuous_alpha, soil_continuous_alpha_1,faith_pd_soil)
soil_continuous_alpha_1 <- data.frame(soil_continuous_alpha_1)
soil_continuous_alpha_1


#for Lime_Buffer_Capacity
glm_shannon_soil_Lime_Buffer_Capacity = glm(Shannon ~ Lime_Buffer_Capacity, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_Lime_Buffer_Capacity)

glm_pd_soil_Lime_Buffer_Capacity = glm(PD ~ Lime_Buffer_Capacity, data=soil_continuous_alpha_1)
summary(glm_pd_soil_Lime_Buffer_Capacity)

#for pH
glm_shannon_soil_pH = glm(Shannon ~ pH, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_pH)

glm_pd_soil_pH = glm(PD ~ pH, data=soil_continuous_alpha_1)
summary(glm_pd_soil_pH)


#for Water_pH
glm_shannon_soil_Water_pH = glm(Shannon ~ Water_pH, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_Water_pH)

glm_pd_soil_Water_pH = glm(PD ~ Water_pH, data=soil_continuous_alpha_1)
summary(glm_pd_soil_Water_pH)


#for Base_Saturation_perc
glm_shannon_soil_Base_Saturation_perc = glm(Shannon ~ Base_Saturation_perc, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_Base_Saturation_perc)

glm_pd_soil_Base_Saturation_perc = glm(PD ~ Base_Saturation_perc, data=soil_continuous_alpha_1)
summary(glm_pd_soil_Base_Saturation_perc)

#for Cation_Exchange_Capacity
glm_shannon_soil_Cation_Exchange_Capacity = glm(Shannon ~ Cation_Exchange_Capacity, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_Cation_Exchange_Capacity)

glm_pd_soil_Cation_Exchange_Capacity = glm(PD ~ Cation_Exchange_Capacity, data=soil_continuous_alpha_1)
summary(glm_pd_soil_Cation_Exchange_Capacity)

#for Ca_ppm
glm_shannon_soil_Ca_ppm = glm(Shannon ~ Ca_ppm, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_Ca_ppm)

glm_pd_soil_Ca_ppm = glm(PD ~ Ca_ppm, data=soil_continuous_alpha_1)
summary(glm_pd_soil_Ca_ppm)

#for Cd_ppm
glm_shannon_soil_Cd_ppm = glm(Shannon ~ Cd_ppm, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_Cd_ppm)

glm_pd_soil_Cd_ppm = glm(PD ~ Cd_ppm, data=soil_continuous_alpha_1)
summary(glm_pd_soil_Cd_ppm)

#for Cr_ppm
glm_shannon_soil_Cr_ppm = glm(Shannon ~ Cr_ppm, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_Cr_ppm)

glm_pd_soil_Cr_ppm = glm(PD ~ Cr_ppm, data=soil_continuous_alpha_1)
summary(glm_pd_soil_Cr_ppm)


#for Cu_ppm
glm_shannon_soil_Cu_ppm = glm(Shannon ~ Cu_ppm, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_Cu_ppm)

glm_pd_soil_Cu_ppm = glm(PD ~ Cu_ppm, data=soil_continuous_alpha_1)
summary(glm_pd_soil_Cu_ppm)

#for Fe_ppm
glm_shannon_soil_Fe_ppm = glm(Shannon ~ Fe_ppm, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_Fe_ppm)

glm_pd_soil_Fe_ppm = glm(PD ~ Fe_ppm, data=soil_continuous_alpha_1)
summary(glm_pd_soil_Fe_ppm)

#for K_ppm
glm_shannon_soil_K_ppm = glm(Shannon ~ K_ppm, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_K_ppm)

glm_pd_soil_K_ppm = glm(PD ~ K_ppm, data=soil_continuous_alpha_1)
summary(glm_pd_soil_K_ppm)

#for Mg_ppm
glm_shannon_soil_Mg_ppm = glm(Shannon ~ Mg_ppm, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_Mg_ppm)

glm_pd_soil_Mg_ppm = glm(PD ~ Mg_ppm, data=soil_continuous_alpha_1)
summary(glm_pd_soil_Mg_ppm)

#for Mn_ppm
glm_shannon_soil_Mn_ppm = glm(Shannon ~ Mn_ppm, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_Mn_ppm)

glm_pd_soil_Mn_ppm = glm(PD ~ Mn_ppm, data=soil_continuous_alpha_1)
summary(glm_pd_soil_Mn_ppm)

#for Mo_ppm
glm_shannon_soil_Mo_ppm = glm(Shannon ~ Mo_ppm, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_Mo_ppm)

glm_pd_soil_Mo_ppm = glm(PD ~ Mo_ppm, data=soil_continuous_alpha_1)
summary(glm_pd_soil_Mo_ppm)

#for Na_ppm
glm_shannon_soil_Na_ppm = glm(Shannon ~ Na_ppm, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_Na_ppm)

glm_pd_soil_Na_ppm = glm(PD ~ Na_ppm, data=soil_continuous_alpha_1)
summary(glm_pd_soil_Na_ppm)

#for Ni_ppm
glm_shannon_soil_Ni_ppm = glm(Shannon ~ Ni_ppm, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_Ni_ppm)

glm_pd_soil_Ni_ppm = glm(PD ~ Ni_ppm, data=soil_continuous_alpha_1)
summary(glm_pd_soil_Ni_ppm)

#for P_ppm
glm_shannon_soil_P_ppm = glm(Shannon ~ P_ppm, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_P_ppm)

glm_pd_soil_P_ppm = glm(PD ~ P_ppm, data=soil_continuous_alpha_1)
summary(glm_pd_soil_P_ppm)

#for Pb_ppm
glm_shannon_soil_Pb_ppm = glm(Shannon ~ Pb_ppm, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_Pb_ppm)

glm_pd_soil_Pb_ppm = glm(PD ~ Pb_ppm, data=soil_continuous_alpha_1)
summary(glm_pd_soil_Pb_ppm)

#for Zn_ppm
glm_shannon_soil_Zn_ppm = glm(Shannon ~ Zn_ppm, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_Zn_ppm)

glm_pd_soil_Zn_ppm = glm(PD ~ Zn_ppm, data=soil_continuous_alpha_1)
summary(glm_pd_soil_Zn_ppm)

#for Organic_Matter_perc
glm_shannon_soil_Organic_Matter_perc = glm(Shannon ~ Organic_Matter_perc, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_Organic_Matter_perc)

glm_pd_soil_Organic_Matter_perc = glm(PD ~ Organic_Matter_perc, data=soil_continuous_alpha_1)
summary(glm_pd_soil_Organic_Matter_perc)

#for N_perc
glm_shannon_soil_N_perc = glm(Shannon ~ N_perc, data=soil_continuous_alpha_1)
summary(glm_shannon_soil_N_perc)

glm_pd_soil_N_perc = glm(PD ~ N_perc, data=soil_continuous_alpha_1)
summary(glm_pd_soil_N_perc)




#For Rhizosphere sample

rhizo_continuous_alpha <- data.frame(bac.css.non_norm_rhizo_rare@sam_data)
rhizo_continuous_alpha
rhizo_continuous_alpha_1 <- estimate_richness(bac.css.non_norm_rhizo_rare, measures = c("Observed", "Shannon"))
rhizo_continuous_alpha_1

# Calculate Faith's Phylogenetic Diversity (PD)
pd_rhizo_otu <- as.data.frame(bac.css.non_norm_rhizo_rare@otu_table)

pd_rhizo_tree <- bac.css.non_norm_rhizo_rare@phy_tree

#check if the tree is rooted
bac.css.non_norm_rhizo_rare@phy_tree

# it is a rooted tree
faith_pd_rhizo <- pd(t(pd_rhizo_otu), pd_rhizo_tree,include.root=T)
faith_pd_rhizo
rhizo_continuous_alpha_1 <- cbind(rhizo_continuous_alpha, rhizo_continuous_alpha_1,faith_pd_rhizo)
rhizo_continuous_alpha_1 <- data.frame(rhizo_continuous_alpha_1)
rhizo_continuous_alpha_1


#for Lime_Buffer_Capacity
glm_shannon_rhizo_Lime_Buffer_Capacity = glm(Shannon ~ Lime_Buffer_Capacity, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_Lime_Buffer_Capacity)

glm_pd_rhizo_Lime_Buffer_Capacity = glm(PD ~ Lime_Buffer_Capacity, data=rhizo_continuous_alpha_1)
summary(glm_pd_rhizo_Lime_Buffer_Capacity)

#for pH
glm_shannon_rhizo_pH = glm(Shannon ~ pH, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_pH)

glm_pd_rhizo_pH = glm(PD ~ pH, data=rhizo_continuous_alpha_1)
summary(glm_pd_rhizo_pH)


#for Water_pH
glm_shannon_rhizo_Water_pH = glm(Shannon ~ Water_pH, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_Water_pH)

glm_pd_rhizo_Water_pH = glm(PD ~ Water_pH, data=rhizo_continuous_alpha_1)
summary(glm_pd_rhizo_Water_pH)


#for Base_Saturation_perc
glm_shannon_rhizo_Base_Saturation_perc = glm(Shannon ~ Base_Saturation_perc, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_Base_Saturation_perc)

glm_pd_rhizo_Base_Saturation_perc = glm(PD ~ Base_Saturation_perc, data=rhizo_continuous_alpha_1)
summary(glm_pd_rhizo_Base_Saturation_perc)

#for Cation_Exchange_Capacity
glm_shannon_rhizo_Cation_Exchange_Capacity = glm(Shannon ~ Cation_Exchange_Capacity, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_Cation_Exchange_Capacity)

glm_pd_rhizo_Cation_Exchange_Capacity = glm(PD ~ Cation_Exchange_Capacity, data=rhizo_continuous_alpha_1)
summary(glm_pd_rhizo_Cation_Exchange_Capacity)

#for Ca_ppm
glm_shannon_rhizo_Ca_ppm = glm(Shannon ~ Ca_ppm, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_Ca_ppm)

glm_pd_rhizo_Ca_ppm = glm(PD ~ Ca_ppm, data=rhizo_continuous_alpha_1)
summary(glm_pd_rhizo_Ca_ppm)

#for Cd_ppm
glm_shannon_rhizo_Cd_ppm = glm(Shannon ~ Cd_ppm, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_Cd_ppm)

glm_pd_rhizo_Cd_ppm = glm(PD ~ Cd_ppm, data=rhizo_continuous_alpha_1)
summary(glm_pd_rhizo_Cd_ppm)

#for Cr_ppm
glm_shannon_rhizo_Cr_ppm = glm(Shannon ~ Cr_ppm, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_Cr_ppm)

glm_pd_rhizo_Cr_ppm = glm(PD ~ Cr_ppm, data=rhizo_continuous_alpha_1)
summary(glm_pd_rhizo_Cr_ppm)


#for Cu_ppm
glm_shannon_rhizo_Cu_ppm = glm(Shannon ~ Cu_ppm, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_Cu_ppm)

glm_pd_rhizo_Cu_ppm = glm(PD ~ Cu_ppm, data=rhizo_continuous_alpha_1)
summary(glm_pd_rhizo_Cu_ppm)

#for Fe_ppm
glm_shannon_rhizo_Fe_ppm = glm(Shannon ~ Fe_ppm, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_Fe_ppm)

glm_pd_rhizo_Fe_ppm = glm(PD ~ Fe_ppm, data=rhizo_continuous_alpha_1)
summary(glm_pd_rhizo_Fe_ppm)

#for K_ppm
glm_shannon_rhizo_K_ppm = glm(Shannon ~ K_ppm, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_K_ppm)

glm_pd_rhizo_K_ppm = glm(PD ~ K_ppm, data=rhizo_continuous_alpha_1)
summary(glm_pd_rhizo_K_ppm)

#for Mg_ppm
glm_shannon_rhizo_Mg_ppm = glm(Shannon ~ Mg_ppm, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_Mg_ppm)

glm_pd_rhizo_Mg_ppm = glm(PD ~ Mg_ppm, data=rhizo_continuous_alpha_1)
summary(glm_pd_rhizo_Mg_ppm)

#for Mn_ppm
glm_shannon_rhizo_Mn_ppm = glm(Shannon ~ Mn_ppm, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_Mn_ppm)

glm_pd_rhizo_Mn_ppm = glm(PD ~ Mn_ppm, data=rhizo_continuous_alpha_1)
summary(glm_pd_rhizo_Mn_ppm)

#for Mo_ppm
glm_shannon_rhizo_Mo_ppm = glm(Shannon ~ Mo_ppm, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_Mo_ppm)

glm_pd_rhizo_Mo_ppm = glm(PD ~ Mo_ppm, data=rhizo_continuous_alpha_1)
summary(glm_pd_rhizo_Mo_ppm)

#for Na_ppm
glm_shannon_rhizo_Na_ppm = glm(Shannon ~ Na_ppm, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_Na_ppm)

glm_pd_rhizo_Na_ppm = glm(PD ~ Na_ppm, data=rhizo_continuous_alpha_1)
summary(glm_pd_rhizo_Na_ppm)

#for Ni_ppm
glm_shannon_rhizo_Ni_ppm = glm(Shannon ~ Ni_ppm, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_Ni_ppm)

glm_pd_rhizo_Ni_ppm = glm(PD ~ Ni_ppm, data=rhizo_continuous_alpha_1)
summary(glm_pd_rhizo_Ni_ppm)

#for P_ppm
glm_shannon_rhizo_P_ppm = glm(Shannon ~ P_ppm, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_P_ppm)

glm_pd_rhizo_P_ppm = glm(PD ~ P_ppm, data=rhizo_continuous_alpha_1)
summary(glm_pd_rhizo_P_ppm)

#for Pb_ppm
glm_shannon_rhizo_Pb_ppm = glm(Shannon ~ Pb_ppm, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_Pb_ppm)

glm_pd_rhizo_Pb_ppm = glm(PD ~ Pb_ppm, data=rhizo_continuous_alpha_1)
summary(glm_pd_rhizo_Pb_ppm)

#for Zn_ppm
glm_shannon_rhizo_Zn_ppm = glm(Shannon ~ Zn_ppm, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_Zn_ppm)

glm_pd_rhizo_Zn_ppm = glm(PD ~ Zn_ppm, data=rhizo_continuous_alpha_1)
summary(glm_pd_rhizo_Zn_ppm)

#for Organic_Matter_perc
glm_shannon_rhizo_Organic_Matter_perc = glm(Shannon ~ Organic_Matter_perc, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_Organic_Matter_perc)

glm_pd_rhizo_Organic_Matter_perc = glm(PD ~ Organic_Matter_perc, data=rhizo_continuous_alpha_1)
summary(glm_pd_rhizo_Organic_Matter_perc)

#for N_perc
glm_shannon_rhizo_N_perc = glm(Shannon ~ N_perc, data=rhizo_continuous_alpha_1)
summary(glm_shannon_rhizo_N_perc)

glm_pd_rhizo_N_perc = glm(PD ~ N_perc, data=rhizo_continuous_alpha_1)
summary(glm_pd_rhizo_N_perc)





#For Root sample

root_continuous_alpha <- data.frame(bac.css.non_norm_root_rare@sam_data)
root_continuous_alpha
root_continuous_alpha_1 <- estimate_richness(bac.css.non_norm_root_rare, measures = c("Observed", "Shannon"))
root_continuous_alpha_1

# Calculate Faith's Phylogenetic Diversity (PD)
pd_root_otu <- as.data.frame(bac.css.non_norm_root_rare@otu_table)

pd_root_tree <- bac.css.non_norm_root_rare@phy_tree

#check if the tree is rooted
bac.css.non_norm_root_rare@phy_tree

# it is a rooted tree
faith_pd_root <- pd(t(pd_root_otu), pd_root_tree,include.root=T)
faith_pd_root
root_continuous_alpha_1 <- cbind(root_continuous_alpha, root_continuous_alpha_1,faith_pd_root)
root_continuous_alpha_1 <- data.frame(root_continuous_alpha_1)
root_continuous_alpha_1



#for Lime_Buffer_Capacity
glm_shannon_root_Lime_Buffer_Capacity = glm(Shannon ~ Lime_Buffer_Capacity, data=root_continuous_alpha_1)
summary(glm_shannon_root_Lime_Buffer_Capacity)

glm_pd_root_Lime_Buffer_Capacity = glm(PD ~ Lime_Buffer_Capacity, data=root_continuous_alpha_1)
summary(glm_pd_root_Lime_Buffer_Capacity)

#for pH
glm_shannon_root_pH = glm(Shannon ~ pH, data=root_continuous_alpha_1)
summary(glm_shannon_root_pH)

glm_pd_root_pH = glm(PD ~ pH, data=root_continuous_alpha_1)
summary(glm_pd_root_pH)


#for Water_pH
glm_shannon_root_Water_pH = glm(Shannon ~ Water_pH, data=root_continuous_alpha_1)
summary(glm_shannon_root_Water_pH)

glm_pd_root_Water_pH = glm(PD ~ Water_pH, data=root_continuous_alpha_1)
summary(glm_pd_root_Water_pH)


#for Base_Saturation_perc
glm_shannon_root_Base_Saturation_perc = glm(Shannon ~ Base_Saturation_perc, data=root_continuous_alpha_1)
summary(glm_shannon_root_Base_Saturation_perc)

glm_pd_root_Base_Saturation_perc = glm(PD ~ Base_Saturation_perc, data=root_continuous_alpha_1)
summary(glm_pd_root_Base_Saturation_perc)

#for Cation_Exchange_Capacity
glm_shannon_root_Cation_Exchange_Capacity = glm(Shannon ~ Cation_Exchange_Capacity, data=root_continuous_alpha_1)
summary(glm_shannon_root_Cation_Exchange_Capacity)

glm_pd_root_Cation_Exchange_Capacity = glm(PD ~ Cation_Exchange_Capacity, data=root_continuous_alpha_1)
summary(glm_pd_root_Cation_Exchange_Capacity)

#for Ca_ppm
glm_shannon_root_Ca_ppm = glm(Shannon ~ Ca_ppm, data=root_continuous_alpha_1)
summary(glm_shannon_root_Ca_ppm)

glm_pd_root_Ca_ppm = glm(PD ~ Ca_ppm, data=root_continuous_alpha_1)
summary(glm_pd_root_Ca_ppm)

#for Cd_ppm
glm_shannon_root_Cd_ppm = glm(Shannon ~ Cd_ppm, data=root_continuous_alpha_1)
summary(glm_shannon_root_Cd_ppm)

glm_pd_root_Cd_ppm = glm(PD ~ Cd_ppm, data=root_continuous_alpha_1)
summary(glm_pd_root_Cd_ppm)

#for Cr_ppm
glm_shannon_root_Cr_ppm = glm(Shannon ~ Cr_ppm, data=root_continuous_alpha_1)
summary(glm_shannon_root_Cr_ppm)

glm_pd_root_Cr_ppm = glm(PD ~ Cr_ppm, data=root_continuous_alpha_1)
summary(glm_pd_root_Cr_ppm)


#for Cu_ppm
glm_shannon_root_Cu_ppm = glm(Shannon ~ Cu_ppm, data=root_continuous_alpha_1)
summary(glm_shannon_root_Cu_ppm)

glm_pd_root_Cu_ppm = glm(PD ~ Cu_ppm, data=root_continuous_alpha_1)
summary(glm_pd_root_Cu_ppm)

#for Fe_ppm
glm_shannon_root_Fe_ppm = glm(Shannon ~ Fe_ppm, data=root_continuous_alpha_1)
summary(glm_shannon_root_Fe_ppm)

glm_pd_root_Fe_ppm = glm(PD ~ Fe_ppm, data=root_continuous_alpha_1)
summary(glm_pd_root_Fe_ppm)

#for K_ppm
glm_shannon_root_K_ppm = glm(Shannon ~ K_ppm, data=root_continuous_alpha_1)
summary(glm_shannon_root_K_ppm)

glm_pd_root_K_ppm = glm(PD ~ K_ppm, data=root_continuous_alpha_1)
summary(glm_pd_root_K_ppm)

#for Mg_ppm
glm_shannon_root_Mg_ppm = glm(Shannon ~ Mg_ppm, data=root_continuous_alpha_1)
summary(glm_shannon_root_Mg_ppm)

glm_pd_root_Mg_ppm = glm(PD ~ Mg_ppm, data=root_continuous_alpha_1)
summary(glm_pd_root_Mg_ppm)

#for Mn_ppm
glm_shannon_root_Mn_ppm = glm(Shannon ~ Mn_ppm, data=root_continuous_alpha_1)
summary(glm_shannon_root_Mn_ppm)

glm_pd_root_Mn_ppm = glm(PD ~ Mn_ppm, data=root_continuous_alpha_1)
summary(glm_pd_root_Mn_ppm)

#for Mo_ppm
glm_shannon_root_Mo_ppm = glm(Shannon ~ Mo_ppm, data=root_continuous_alpha_1)
summary(glm_shannon_root_Mo_ppm)

glm_pd_root_Mo_ppm = glm(PD ~ Mo_ppm, data=root_continuous_alpha_1)
summary(glm_pd_root_Mo_ppm)

#for Na_ppm
glm_shannon_root_Na_ppm = glm(Shannon ~ Na_ppm, data=root_continuous_alpha_1)
summary(glm_shannon_root_Na_ppm)

glm_pd_root_Na_ppm = glm(PD ~ Na_ppm, data=root_continuous_alpha_1)
summary(glm_pd_root_Na_ppm)

#for Ni_ppm
glm_shannon_root_Ni_ppm = glm(Shannon ~ Ni_ppm, data=root_continuous_alpha_1)
summary(glm_shannon_root_Ni_ppm)

glm_pd_root_Ni_ppm = glm(PD ~ Ni_ppm, data=root_continuous_alpha_1)
summary(glm_pd_root_Ni_ppm)

#for P_ppm
glm_shannon_root_P_ppm = glm(Shannon ~ P_ppm, data=root_continuous_alpha_1)
summary(glm_shannon_root_P_ppm)

glm_pd_root_P_ppm = glm(PD ~ P_ppm, data=root_continuous_alpha_1)
summary(glm_pd_root_P_ppm)

#for Pb_ppm
glm_shannon_root_Pb_ppm = glm(Shannon ~ Pb_ppm, data=root_continuous_alpha_1)
summary(glm_shannon_root_Pb_ppm)

glm_pd_root_Pb_ppm = glm(PD ~ Pb_ppm, data=root_continuous_alpha_1)
summary(glm_pd_root_Pb_ppm)

#for Zn_ppm
glm_shannon_root_Zn_ppm = glm(Shannon ~ Zn_ppm, data=root_continuous_alpha_1)
summary(glm_shannon_root_Zn_ppm)

glm_pd_root_Zn_ppm = glm(PD ~ Zn_ppm, data=root_continuous_alpha_1)
summary(glm_pd_root_Zn_ppm)

#for Organic_Matter_perc
glm_shannon_root_Organic_Matter_perc = glm(Shannon ~ Organic_Matter_perc, data=root_continuous_alpha_1)
summary(glm_shannon_root_Organic_Matter_perc)

glm_pd_root_Organic_Matter_perc = glm(PD ~ Organic_Matter_perc, data=root_continuous_alpha_1)
summary(glm_pd_root_Organic_Matter_perc)

#for N_perc
glm_shannon_root_N_perc = glm(Shannon ~ N_perc, data=root_continuous_alpha_1)
summary(glm_shannon_root_N_perc)

glm_pd_root_N_perc = glm(PD ~ N_perc, data=root_continuous_alpha_1)
summary(glm_pd_root_N_perc)




#For Stalk sample

stalk_continuous_alpha <- data.frame(bac.css.non_norm_stalk_rare@sam_data)
stalk_continuous_alpha
stalk_continuous_alpha_1 <- estimate_richness(bac.css.non_norm_stalk_rare, measures = c("Observed", "Shannon"))
stalk_continuous_alpha_1

# Calculate Faith's Phylogenetic Diversity (PD)
pd_stalk_otu <- as.data.frame(bac.css.non_norm_stalk_rare@otu_table)

pd_stalk_tree <- bac.css.non_norm_stalk_rare@phy_tree

#check if the tree is rooted
bac.css.non_norm_stalk_rare@phy_tree

# it is a rooted tree
faith_pd_stalk <- pd(t(pd_stalk_otu), pd_stalk_tree,include.root=T)
faith_pd_stalk
stalk_continuous_alpha_1 <- cbind(stalk_continuous_alpha, stalk_continuous_alpha_1,faith_pd_stalk)
stalk_continuous_alpha_1 <- data.frame(stalk_continuous_alpha_1)
stalk_continuous_alpha_1





#for Lime_Buffer_Capacity
glm_shannon_stalk_Lime_Buffer_Capacity = glm(Shannon ~ Lime_Buffer_Capacity, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_Lime_Buffer_Capacity)

glm_pd_stalk_Lime_Buffer_Capacity = glm(PD ~ Lime_Buffer_Capacity, data=stalk_continuous_alpha_1)
summary(glm_pd_stalk_Lime_Buffer_Capacity)

#for pH
glm_shannon_stalk_pH = glm(Shannon ~ pH, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_pH)

glm_pd_stalk_pH = glm(PD ~ pH, data=stalk_continuous_alpha_1)
summary(glm_pd_stalk_pH)


#for Water_pH
glm_shannon_stalk_Water_pH = glm(Shannon ~ Water_pH, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_Water_pH)

glm_pd_stalk_Water_pH = glm(PD ~ Water_pH, data=stalk_continuous_alpha_1)
summary(glm_pd_stalk_Water_pH)


#for Base_Saturation_perc
glm_shannon_stalk_Base_Saturation_perc = glm(Shannon ~ Base_Saturation_perc, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_Base_Saturation_perc)

glm_pd_stalk_Base_Saturation_perc = glm(PD ~ Base_Saturation_perc, data=stalk_continuous_alpha_1)
summary(glm_pd_stalk_Base_Saturation_perc)

#for Cation_Exchange_Capacity
glm_shannon_stalk_Cation_Exchange_Capacity = glm(Shannon ~ Cation_Exchange_Capacity, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_Cation_Exchange_Capacity)

glm_pd_stalk_Cation_Exchange_Capacity = glm(PD ~ Cation_Exchange_Capacity, data=stalk_continuous_alpha_1)
summary(glm_pd_stalk_Cation_Exchange_Capacity)

#for Ca_ppm
glm_shannon_stalk_Ca_ppm = glm(Shannon ~ Ca_ppm, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_Ca_ppm)

glm_pd_stalk_Ca_ppm = glm(PD ~ Ca_ppm, data=stalk_continuous_alpha_1)
summary(glm_pd_stalk_Ca_ppm)

#for Cd_ppm
glm_shannon_stalk_Cd_ppm = glm(Shannon ~ Cd_ppm, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_Cd_ppm)

glm_pd_stalk_Cd_ppm = glm(PD ~ Cd_ppm, data=stalk_continuous_alpha_1)
summary(glm_pd_stalk_Cd_ppm)

#for Cr_ppm
glm_shannon_stalk_Cr_ppm = glm(Shannon ~ Cr_ppm, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_Cr_ppm)

glm_pd_stalk_Cr_ppm = glm(PD ~ Cr_ppm, data=stalk_continuous_alpha_1)
summary(glm_pd_stalk_Cr_ppm)


#for Cu_ppm
glm_shannon_stalk_Cu_ppm = glm(Shannon ~ Cu_ppm, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_Cu_ppm)

glm_pd_stalk_Cu_ppm = glm(PD ~ Cu_ppm, data=stalk_continuous_alpha_1)
summary(glm_pd_stalk_Cu_ppm)

#for Fe_ppm
glm_shannon_stalk_Fe_ppm = glm(Shannon ~ Fe_ppm, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_Fe_ppm)

glm_pd_stalk_Fe_ppm = glm(PD ~ Fe_ppm, data=stalk_continuous_alpha_1)
summary(glm_pd_stalk_Fe_ppm)

#for K_ppm
glm_shannon_stalk_K_ppm = glm(Shannon ~ K_ppm, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_K_ppm)

glm_pd_stalk_K_ppm = glm(PD ~ K_ppm, data=stalk_continuous_alpha_1)
summary(glm_pd_stalk_K_ppm)

#for Mg_ppm
glm_shannon_stalk_Mg_ppm = glm(Shannon ~ Mg_ppm, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_Mg_ppm)

glm_pd_stalk_Mg_ppm = glm(PD ~ Mg_ppm, data=stalk_continuous_alpha_1)
summary(glm_pd_stalk_Mg_ppm)

#for Mn_ppm
glm_shannon_stalk_Mn_ppm = glm(Shannon ~ Mn_ppm, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_Mn_ppm)

glm_pd_stalk_Mn_ppm = glm(PD ~ Mn_ppm, data=stalk_continuous_alpha_1)
summary(glm_pd_stalk_Mn_ppm)

#for Mo_ppm
glm_shannon_stalk_Mo_ppm = glm(Shannon ~ Mo_ppm, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_Mo_ppm)

glm_pd_stalk_Mo_ppm = glm(PD ~ Mo_ppm, data=stalk_continuous_alpha_1)
summary(glm_pd_stalk_Mo_ppm)

#for Na_ppm
glm_shannon_stalk_Na_ppm = glm(Shannon ~ Na_ppm, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_Na_ppm)

glm_pd_stalk_Na_ppm = glm(PD ~ Na_ppm, data=stalk_continuous_alpha_1)
summary(glm_pd_stalk_Na_ppm)

#for Ni_ppm
glm_shannon_stalk_Ni_ppm = glm(Shannon ~ Ni_ppm, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_Ni_ppm)

glm_pd_stalk_Ni_ppm = glm(PD ~ Ni_ppm, data=stalk_continuous_alpha_1)
summary(glm_pd_stalk_Ni_ppm)

#for P_ppm
glm_shannon_stalk_P_ppm = glm(Shannon ~ P_ppm, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_P_ppm)

glm_pd_stalk_P_ppm = glm(PD ~ P_ppm, data=stalk_continuous_alpha_1)
summary(glm_pd_stalk_P_ppm)

#for Pb_ppm
glm_shannon_stalk_Pb_ppm = glm(Shannon ~ Pb_ppm, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_Pb_ppm)

glm_pd_stalk_Pb_ppm = glm(PD ~ Pb_ppm, data=stalk_continuous_alpha_1)
summary(glm_pd_stalk_Pb_ppm)

#for Zn_ppm
glm_shannon_stalk_Zn_ppm = glm(Shannon ~ Zn_ppm, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_Zn_ppm)

glm_pd_stalk_Zn_ppm = glm(PD ~ Zn_ppm, data=stalk_continuous_alpha_1)
summary(glm_pd_stalk_Zn_ppm)

#for Organic_Matter_perc
glm_shannon_stalk_Organic_Matter_perc = glm(Shannon ~ Organic_Matter_perc, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_Organic_Matter_perc)

glm_pd_stalk_Organic_Matter_perc = glm(PD ~ Organic_Matter_perc, data=stalk_continuous_alpha_1)
summary(glm_pd_stalk_Organic_Matter_perc)

#for N_perc
glm_shannon_stalk_N_perc = glm(Shannon ~ N_perc, data=stalk_continuous_alpha_1)
summary(glm_shannon_stalk_N_perc)

glm_pd_stalk_N_perc = glm(PD ~ N_perc, data=stalk_continuous_alpha_1)
summary(glm_pd_stalk_N_perc)
