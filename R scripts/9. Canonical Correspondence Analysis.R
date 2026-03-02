#Load required libraries
library(ggplot2)
library(phyloseq)
library(dplyr)
library(vegan)
library(car)  # For VIF analysis


set.seed(1234)

# ================================
# read in the non-normalized data
# ================================
rhizo <- readRDS("C:/Aduragbemi/Microbiome/RDS/rhizo_bac_sub_no_bad_Filtered_final_non_normalized.rds")
soil <- readRDS("C:/Aduragbemi/Microbiome/RDS/soil_bac_sub_no_bad_Filtered_final_non_normalized.rds")
root <- readRDS("C:/Aduragbemi/Microbiome/RDS/root_bac_sub_no_bad_Filtered_final_non_normalized.rds")



#For Soil Sample
# Convert sample metadata from phyloseq class to data frame
sample_data_soil_df = data.frame(sample_data(soil))

# Extract geographic information
geo_soil <- sample_data_soil_df[, c("GPSlatitude", "GPSlongitude")]

# Extract environmental variables
env_soil <- sample_data_soil_df[, c("Growing_Degree_Days", "precipIntensity", "precipProbability", "temperatureHigh",
                                    "temperatureLow", "dewPoint", "humidity", "windSpeed", "uvIndex")]

# Extract soil properties
prop_soil <- sample_data_soil_df[, c("Lime_Buffer_Capacity", "pH", "Water_pH", "Base_Saturation_perc",
                                     "Cation_Exchange_Capacity", "Ca_ppm", "Cu_ppm", "Fe_ppm", "K_ppm",
                                     "Mg_ppm", "Mn_ppm", "Na_ppm", "Ni_ppm", "P_ppm","Pb_ppm", "Zn_ppm", 
                                     "Organic_Matter_perc", "N_perc", "Mo_ppm", "Cd_ppm", "Cr_ppm")]

# Check for NA values
print("Are there NA values in soil properties:")
print(any(is.na(prop_soil)))

# Rename column names
env_soil <- env_soil %>%
  rename(
    "Growing Degree Days" = "Growing_Degree_Days", "PI" = "precipIntensity", "PP" ="precipProbability", "HT" = "temperatureHigh",
    "LT" = "temperatureLow", "DP" = "dewPoint", "HUM" = "humidity", "WS" = "windSpeed", "UVI" = "uvIndex")

geo_soil <- geo_soil %>%
  rename("LAT" = "GPSlatitude", "LONG" = "GPSlongitude")

prop_soil <- prop_soil %>%
  rename(
    "LBC" = "Lime_Buffer_Capacity", "pH" = "pH", "WpH" = "Water_pH", "BS" = "Base_Saturation_perc",
    "CEC" = "Cation_Exchange_Capacity", "Ca" = "Ca_ppm", "Cu" = "Cu_ppm", "Fe" = "Fe_ppm", "K" = "K_ppm",
    "Mg" = "Mg_ppm", "Mn" = "Mn_ppm", "Na" = "Na_ppm", "Ni" = "Ni_ppm", "P" = "P_ppm","Pb" = "Pb_ppm", "Zn" = "Zn_ppm", 
    "OM" = "Organic_Matter_perc", "N" = "N_perc", "Mo" = "Mo_ppm", "Cd" = "Cd_ppm", "Cr" = "Cr_ppm")


# Get ASV table
asv_soil <- as.data.frame(otu_table(soil))
asv_soil_t <- t(asv_soil)

#Hellinger transformation on species abundance data 
#This turns absolute abundance into relative abundance
#We will call this spec.h
asv_soil_t <- decostand(asv_soil_t, method = "hellinger")

#check whether to use CCA or RDA
decorana(asv_soil_t)

#the length of the first axis (DCA1) is 13.2294, which is greater than 3. 
#This suggests that CCA would be more appropriate


# Check for multicollinearity in environmental variables
# Create a data frame with all environmental variables
env_model_soil <- lm(`Growing Degree Days` ~ ., data = env_soil)
vif_env_soil <- vif(env_model_soil)

# Identify variables with VIF > 10
high_vif_env_soil <- names(vif_env_soil[vif_env_soil > 10])

# Check for multicollinearity in soil properties
# Due to many variables, use stepwise regression to avoid perfect multicollinearity
prop_model_soil <- lm(pH ~ ., data = prop_soil)
vif_prop_soil <- vif(prop_model_soil)
vif_prop_soil

# Identify variables with VIF > 10
high_vif_prop_soil <- names(vif_prop_soil[vif_prop_soil > 10])

#Let's get rid of WpH, LT, DP, pH, Ca, , CEC, Mg, PP, OMP 
# Remove high-VIF variables from soil properties data frame
prop_soil_vif_filtered <- prop_soil[, 
                                      !colnames(prop_soil) %in% high_vif_prop_soil, 
                                      drop = FALSE]

# Check remaining variables
colnames(prop_soil_vif_filtered)

# Method 1: CCA analysis using only environmental variables after removing collinear variables
spec.cca_env_soil <- cca(asv_soil_t ~ ., data = env_soil)
print(summary(spec.cca_env_soil))

# Method 2: CCA analysis using only soil properties
spec.cca_prop_soil <- cca(asv_soil_t ~ ., data = prop_soil_vif_filtered)
print(summary(spec.cca_prop_soil))

# Method 3: CCA analysis using only geographic information
spec.cca_geo_soil <- cca(asv_soil_t ~ ., data = geo_soil)
print(summary(spec.cca_geo_soil))

# Method 4: Variance partitioning analysis to compare contributions of three variable groups
vp_soil <- varpart(asv_soil_t, env_soil, prop_soil, geo_soil)
vp_soil

# Plot variance partitioning diagram
plot(vp_soil, Xnames = c("Environment", "Soil Properties", "Geography"), bg = c("red", "blue", "green"))


# Use forward selection to select important variables from environmental variables
step.env_soil <- ordistep(cca(asv_soil_t ~ 1, data = env_soil), 
                     scope = formula(spec.cca_env_soil), 
                     direction = "forward", 
                     permutations = 999)
step.env_soil 
# Use forward selection to select important variables from soil properties
step.prop_soil <- ordistep(cca(asv_soil_t ~ 1, data = prop_soil), 
                      scope = formula(spec.cca_prop_soil), 
                      direction = "forward", 
                      permutations = 999)
step.prop_soil
# Use forward selection to select important variables from geographic information
step.geo_soil <- ordistep(cca(asv_soil_t ~ 1, data = geo_soil), 
                     scope = formula(spec.cca_geo_soil), 
                     direction = "forward", 
                     permutations = 999)
step.geo_soil
# Extract selected variables
selected_env_vars_soil <- attr(terms(step.env_soil$call$formula), "term.labels")
selected_prop_vars_soil <- attr(terms(step.prop_soil$call$formula), "term.labels")
selected_geo_vars_soil <- attr(terms(step.geo_soil$call$formula), "term.labels")

# Create a data frame with all selected variables
selected_vars_soil <- data.frame(
  env_soil[, intersect(selected_env_vars_soil, colnames(env_soil)), drop = FALSE],
  prop_soil[, intersect(selected_prop_vars_soil, colnames(prop_soil)), drop = FALSE],
  geo_soil[, intersect(selected_geo_vars_soil, colnames(geo_soil)), drop = FALSE]
)
selected_vars_soil
# Perform integrated CCA analysis with all selected variables
spec.cca_combined_soil <- cca(asv_soil_t ~ ., data = selected_vars_soil)
print(summary(spec.cca_combined_soil))

# Extract CCA results
cca_combined_soil <- scores(spec.cca_combined_soil, display = c("sites", "species", "bp"), scaling = 1)

# Extract site and environmental variable coordinates
sites_combined_soil <- as.data.frame(cca_combined_soil$sites)
bp_combined_soil <- as.data.frame(cca_combined_soil$biplot)

# Merge site coordinates with metadata
samples_combined_soil <- cbind(sites_combined_soil, sample_data_soil_df)

# Calculate percentage of variance explained by CCA axes
aa_combined_soil <- summary(spec.cca_combined_soil)$cont$importance[2, 1:2]
aa_combined_soil <- rbind(aa_combined_soil, aa_combined_soil * 100)

# Create environmental variables data frame for plotting
envi_combined_soil <- data.frame(
  env = rownames(bp_combined_soil),
  CCA1 = bp_combined_soil$CCA1,
  CCA2 = bp_combined_soil$CCA2
)

# Add variable type colors
envi_combined_soil$type <- "Other"
envi_combined_soil$type[envi_combined_soil$env %in% selected_env_vars_soil] <- "Environment"
envi_combined_soil$type[envi_combined_soil$env %in% selected_prop_vars_soil] <- "soil"
envi_combined_soil$type[envi_combined_soil$env %in% selected_geo_vars_soil] <- "Geography"

# Plot integrated CCA diagram
cca_soil <- ggplot() + 
  geom_segment(data = envi_combined_soil, aes(x = 0, y = 0, xend = CCA1*3.1, yend = CCA2*3.1, color = type), 
               arrow = arrow(length = unit(0.5, "cm")), linewidth = 1.2) + 
  geom_text(data = envi_combined_soil, aes(x = CCA1*3.1, y = CCA2*3.1, label = env, color = type), 
            size = 10, check_overlap = FALSE) +
  geom_point(data = samples_combined_soil, aes(x=CCA1, y=CCA2, shape = factor(`Growing_Degree_Days`)), size = 7) +
  scale_color_manual(values = c("Environment" = "red", "soil" = "blue", "Geography" = "green", "Other" = "gray")) +
  theme_bw() + 
  theme(panel.grid.major=element_line(colour=NA), panel.grid.minor = element_blank()) + 
  xlab(paste('CCA1 (', aa_combined_soil[2,1], '%)', sep ='')) + 
  ylab(paste('CCA2 (', aa_combined_soil[2,2], '%)', sep ='')) + 
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", face = "bold", size = 15),      
        axis.text.x=element_text(angle=0,hjust=1,vjust=1.0,size=15, colour = "black"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black"),
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 10, face = "bold", angle = 0),
        legend.position= "right",
        legend.title = element_text(face="bold", size = 15)) +
  labs(color = "Variable Type", shape = "Growing Degree Days") +
  ggtitle("Soil Microbiome CCA Analysis")

cca_soil

#save plot in pdf and png
ggsave("cca_soil.png", plot = cca_soil , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("cca_soil.pdf", plot = cca_soil, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")


# Test integrated model significance
# ANOVA - test significance of the entire CCA model
anova_combined_soil <- anova(spec.cca_combined_soil, perm.max=999)
anova_combined_soil
# Test significance of each variable
anova_terms_soil <- anova(spec.cca_combined_soil, by="terms", perm.max=999)
anova_terms_soil
# Extract raw p-values excluding NA values
raw_pvals_soil <- anova_terms_soil$`Pr(>F)`[!is.na(anova_terms_soil$`Pr(>F)`)]

# Apply Bonferroni correction
bonf_pvals_soil <- p.adjust(raw_pvals_soil, method = "bonferroni")

# Apply other correction methods for comparison
fdr_pvals_soil <- p.adjust(raw_pvals_soil, method = "fdr")
holm_pvals_soil <- p.adjust(raw_pvals_soil, method = "holm")

# Create complete correction results table
correction_results_soil <- data.frame(
  Variable = rownames(anova_terms_soil)[!is.na(anova_terms_soil$`Pr(>F)`)],
  F_value = anova_terms_soil$F[!is.na(anova_terms_soil$`Pr(>F)`)],
  Raw_P_value = raw_pvals_soil,
  Bonferroni_P_value = bonf_pvals_soil,
  FDR_P_value = fdr_pvals_soil,
  Holm_P_value = holm_pvals_soil
)

# Add variable type information
correction_results_soil$Type <- "Other"
correction_results_soil$Type[correction_results_soil$Variable %in% selected_env_vars_soil] <- "Environment"
correction_results_soil$Type[correction_results_soil$Variable %in% selected_prop_vars_soil] <- "soil"
correction_results_soil$Type[correction_results_soil$Variable %in% selected_geo_vars_soil] <- "Geography"

# Calculate adjusted R-squared
rsq_soil <- RsquareAdj(spec.cca_combined_soil)

# Save results table to CSV file
write.csv(correction_results_soil, "C:/Aduragbemi/Manuscript/Review 5/PDF5/cca_correction_results_soil.csv", row.names = FALSE)




#For Rhizosphere Sample
# Convert sample metadata from phyloseq class to data frame
sample_data_rhizo_df = data.frame(sample_data(rhizo))

# Extract geographic information
geo_rhizo <- sample_data_rhizo_df[, c("GPSlatitude", "GPSlongitude")]

# Extract environmental variables
env_rhizo <- sample_data_rhizo_df[, c("Growing_Degree_Days", "precipIntensity", "precipProbability", "temperatureHigh",
                                    "temperatureLow", "dewPoint", "humidity", "windSpeed", "uvIndex")]

# Extract rhizo properties
prop_rhizo <- sample_data_rhizo_df[, c("Lime_Buffer_Capacity", "pH", "Water_pH", "Base_Saturation_perc",
                                     "Cation_Exchange_Capacity", "Ca_ppm", "Cu_ppm", "Fe_ppm", "K_ppm",
                                     "Mg_ppm", "Mn_ppm", "Na_ppm", "Ni_ppm", "P_ppm","Pb_ppm", "Zn_ppm", 
                                     "Organic_Matter_perc", "N_perc", "Mo_ppm", "Cd_ppm", "Cr_ppm")]

# Check for NA values
print("Are there NA values in rhizo properties:")
print(any(is.na(prop_rhizo)))

# Rename column names
env_rhizo <- env_rhizo %>%
  rename(
    "Growing Degree Days" = "Growing_Degree_Days", "PI" = "precipIntensity", "PP" ="precipProbability", "HT" = "temperatureHigh",
    "LT" = "temperatureLow", "DP" = "dewPoint", "HUM" = "humidity", "WS" = "windSpeed", "UVI" = "uvIndex")

geo_rhizo <- geo_rhizo %>%
  rename("LAT" = "GPSlatitude", "LONG" = "GPSlongitude")

prop_rhizo <- prop_rhizo %>%
  rename(
    "LBC" = "Lime_Buffer_Capacity", "pH" = "pH", "WpH" = "Water_pH", "BS" = "Base_Saturation_perc",
    "CEC" = "Cation_Exchange_Capacity", "Ca" = "Ca_ppm", "Cu" = "Cu_ppm", "Fe" = "Fe_ppm", "K" = "K_ppm",
    "Mg" = "Mg_ppm", "Mn" = "Mn_ppm", "Na" = "Na_ppm", "Ni" = "Ni_ppm", "P" = "P_ppm","Pb" = "Pb_ppm", "Zn" = "Zn_ppm", 
    "OM" = "Organic_Matter_perc", "N" = "N_perc", "Mo" = "Mo_ppm", "Cd" = "Cd_ppm", "Cr" = "Cr_ppm")


# Get ASV table
asv_rhizo <- as.data.frame(otu_table(rhizo))
asv_rhizo_t <- t(asv_rhizo)

#Hellinger transformation on species abundance data 
#This turns absolute abundance into relative abundance
#We will call this spec.h
asv_rhizo_t <- decostand(asv_rhizo_t, method = "hellinger")

#check whether to use CCA or RDA
decorana(asv_rhizo_t)

#the length of the first axis (DCA1) is 13.2294, which is greater than 3. 
#This suggests that CCA would be more appropriate


# Check for multicollinearity in environmental variables
# Create a data frame with all environmental variables
env_model_rhizo <- lm(`Growing Degree Days` ~ ., data = env_rhizo)
vif_env_rhizo <- vif(env_model_rhizo)

# Identify variables with VIF > 10
high_vif_env_rhizo <- names(vif_env_rhizo[vif_env_rhizo > 10])

# Check for multicollinearity in rhizo properties
# Due to many variables, use stepwise regression to avoid perfect multicollinearity
prop_model_rhizo <- lm(pH ~ ., data = prop_rhizo)
vif_prop_rhizo <- vif(prop_model_rhizo)
vif_prop_rhizo

# Identify variables with VIF > 10
high_vif_prop_rhizo <- names(vif_prop_rhizo[vif_prop_rhizo > 10])

#Let's get rid of WpH, LT, DP, pH, Ca, , CEC, Mg, PP, OMP 
# Remove high-VIF variables from rhizo properties data frame
prop_rhizo_vif_filtered <- prop_rhizo[, 
                                    !colnames(prop_rhizo) %in% high_vif_prop_rhizo, 
                                    drop = FALSE]

# Check remaining variables
colnames(prop_rhizo_vif_filtered)

# Method 1: CCA analysis using only environmental variables after removing collinear variables
spec.cca_env_rhizo <- cca(asv_rhizo_t ~ ., data = env_rhizo)
print(summary(spec.cca_env_rhizo))

# Method 2: CCA analysis using only rhizo properties
spec.cca_prop_rhizo <- cca(asv_rhizo_t ~ ., data = prop_rhizo_vif_filtered)
print(summary(spec.cca_prop_rhizo))

# Method 3: CCA analysis using only geographic information
spec.cca_geo_rhizo <- cca(asv_rhizo_t ~ ., data = geo_rhizo)
print(summary(spec.cca_geo_rhizo))

# Method 4: Variance partitioning analysis to compare contributions of three variable groups
vp_rhizo <- varpart(asv_rhizo_t, env_rhizo, prop_rhizo, geo_rhizo)
vp_rhizo

# Plot variance partitioning diagram
plot(vp_rhizo, Xnames = c("Environment", "rhizo Properties", "Geography"), bg = c("red", "blue", "green"))


# Use forward selection to select important variables from environmental variables
step.env_rhizo <- ordistep(cca(asv_rhizo_t ~ 1, data = env_rhizo), 
                          scope = formula(spec.cca_env_rhizo), 
                          direction = "forward", 
                          permutations = 999)

# Use forward selection to select important variables from rhizo properties
step.prop_rhizo <- ordistep(cca(asv_rhizo_t ~ 1, data = prop_rhizo), 
                           scope = formula(spec.cca_prop_rhizo), 
                           direction = "forward", 
                           permutations = 999)

# Use forward selection to select important variables from geographic information
step.geo_rhizo <- ordistep(cca(asv_rhizo_t ~ 1, data = geo_rhizo), 
                          scope = formula(spec.cca_geo_rhizo), 
                          direction = "forward", 
                          permutations = 999)

# Extract selected variables
selected_env_vars_rhizo <- attr(terms(step.env_rhizo$call$formula), "term.labels")
selected_prop_vars_rhizo <- attr(terms(step.prop_rhizo$call$formula), "term.labels")
selected_geo_vars_rhizo <- attr(terms(step.geo_rhizo$call$formula), "term.labels")

# Create a data frame with all selected variables
selected_vars_rhizo <- data.frame(
  env_rhizo[, intersect(selected_env_vars_rhizo, colnames(env_rhizo)), drop = FALSE],
  prop_rhizo[, intersect(selected_prop_vars_rhizo, colnames(prop_rhizo)), drop = FALSE],
  geo_rhizo[, intersect(selected_geo_vars_rhizo, colnames(geo_rhizo)), drop = FALSE]
)

# Perform integrated CCA analysis with all selected variables
spec.cca_combined_rhizo <- cca(asv_rhizo_t ~ ., data = selected_vars_rhizo)
print(summary(spec.cca_combined_rhizo))

# Extract CCA results
cca_combined_rhizo <- scores(spec.cca_combined_rhizo, display = c("sites", "species", "bp"), scaling = 1)

# Extract site and environmental variable coordinates
sites_combined_rhizo <- as.data.frame(cca_combined_rhizo$sites)
bp_combined_rhizo <- as.data.frame(cca_combined_rhizo$biplot)

# Merge site coordinates with metadata
samples_combined_rhizo <- cbind(sites_combined_rhizo, sample_data_rhizo_df)

# Calculate percentage of variance explained by CCA axes
aa_combined_rhizo <- summary(spec.cca_combined_rhizo)$cont$importance[2, 1:2]
aa_combined_rhizo <- rbind(aa_combined_rhizo, aa_combined_rhizo * 100)

# Create environmental variables data frame for plotting
envi_combined_rhizo <- data.frame(
  env = rownames(bp_combined_rhizo),
  CCA1 = bp_combined_rhizo$CCA1,
  CCA2 = bp_combined_rhizo$CCA2
)

# Add variable type colors
envi_combined_rhizo$type <- "Other"
envi_combined_rhizo$type[envi_combined_rhizo$env %in% selected_env_vars_rhizo] <- "Environment"
envi_combined_rhizo$type[envi_combined_rhizo$env %in% selected_prop_vars_rhizo] <- "rhizo"
envi_combined_rhizo$type[envi_combined_rhizo$env %in% selected_geo_vars_rhizo] <- "Geography"

# Plot integrated CCA diagram
cca_rhizo <- ggplot() + 
  geom_segment(data = envi_combined_rhizo, aes(x = 0, y = 0, xend = CCA1*3.1, yend = CCA2*3.1, color = type), 
               arrow = arrow(length = unit(0.5, "cm")), linewidth = 1.2) + 
  geom_text(data = envi_combined_rhizo, aes(x = CCA1*3.1, y = CCA2*3.1, label = env, color = type), 
            size = 10, check_overlap = FALSE) +
  geom_point(data = samples_combined_rhizo, aes(x=CCA1, y=CCA2, shape = factor(`Growing_Degree_Days`)), size = 7) +
  scale_color_manual(values = c("Environment" = "red", "rhizo" = "blue", "Geography" = "green", "Other" = "gray")) +
  theme_bw() + 
  theme(panel.grid.major=element_line(colour=NA), panel.grid.minor = element_blank()) + 
  xlab(paste('CCA1 (', aa_combined_rhizo[2,1], '%)', sep ='')) + 
  ylab(paste('CCA2 (', aa_combined_rhizo[2,2], '%)', sep ='')) + 
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", face = "bold", size = 15),      
        axis.text.x=element_text(angle=0,hjust=1,vjust=1.0,size=15, colour = "black"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black"),
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 10, face = "bold", angle = 0),
        legend.position= "right",
        legend.title = element_text(face="bold", size = 15)) +
  labs(color = "Variable Type", shape = "Growing Degree Days") +
  ggtitle("rhizo Microbiome CCA Analysis")

cca_rhizo

#save plot in pdf and png
ggsave("cca_rhizo.png", plot = cca_rhizo , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("cca_rhizo.pdf", plot = cca_rhizo, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")


# Test integrated model significance
# ANOVA - test significance of the entire CCA model
anova_combined_rhizo <- anova(spec.cca_combined_rhizo, perm.max=999)

# Test significance of each variable
anova_terms_rhizo <- anova(spec.cca_combined_rhizo, by="terms", perm.max=999)

# Extract raw p-values excluding NA values
raw_pvals_rhizo <- anova_terms_rhizo$`Pr(>F)`[!is.na(anova_terms_rhizo$`Pr(>F)`)]

# Apply Bonferroni correction
bonf_pvals_rhizo <- p.adjust(raw_pvals_rhizo, method = "bonferroni")

# Apply other correction methods for comparison
fdr_pvals_rhizo <- p.adjust(raw_pvals_rhizo, method = "fdr")
holm_pvals_rhizo <- p.adjust(raw_pvals_rhizo, method = "holm")

# Create complete correction results table
correction_results_rhizo <- data.frame(
  Variable = rownames(anova_terms_rhizo)[!is.na(anova_terms_rhizo$`Pr(>F)`)],
  F_value = anova_terms_rhizo$F[!is.na(anova_terms_rhizo$`Pr(>F)`)],
  Raw_P_value = raw_pvals_rhizo,
  Bonferroni_P_value = bonf_pvals_rhizo,
  FDR_P_value = fdr_pvals_rhizo,
  Holm_P_value = holm_pvals_rhizo
)

# Add variable type information
correction_results_rhizo$Type <- "Other"
correction_results_rhizo$Type[correction_results_rhizo$Variable %in% selected_env_vars_rhizo] <- "Environment"
correction_results_rhizo$Type[correction_results_rhizo$Variable %in% selected_prop_vars_rhizo] <- "rhizo"
correction_results_rhizo$Type[correction_results_rhizo$Variable %in% selected_geo_vars_rhizo] <- "Geography"

# Calculate adjusted R-squared
rsq_rhizo <- RsquareAdj(spec.cca_combined_rhizo)

# Save results table to CSV file
write.csv(correction_results_rhizo, "C:/Aduragbemi/Manuscript/Review 5/PDF5/cca_correction_results_rhizo.csv", row.names = FALSE)





#For Root Sample
# Convert sample metadata from phyloseq class to data frame
sample_data_root_df = data.frame(sample_data(root))

# Extract geographic information
geo_root <- sample_data_root_df[, c("GPSlatitude", "GPSlongitude")]

# Extract environmental variables
env_root <- sample_data_root_df[, c("Growing_Degree_Days", "precipIntensity", "precipProbability", "temperatureHigh",
                                      "temperatureLow", "dewPoint", "humidity", "windSpeed", "uvIndex")]

# Extract root properties
prop_root <- sample_data_root_df[, c("Lime_Buffer_Capacity", "pH", "Water_pH", "Base_Saturation_perc",
                                       "Cation_Exchange_Capacity", "Ca_ppm", "Cu_ppm", "Fe_ppm", "K_ppm",
                                       "Mg_ppm", "Mn_ppm", "Na_ppm", "Ni_ppm", "P_ppm","Pb_ppm", "Zn_ppm", 
                                       "Organic_Matter_perc", "N_perc", "Mo_ppm", "Cd_ppm", "Cr_ppm")]

# Check for NA values
print("Are there NA values in root properties:")
print(any(is.na(prop_root)))

# Rename column names
env_root <- env_root %>%
  rename(
    "Growing Degree Days" = "Growing_Degree_Days", "PI" = "precipIntensity", "PP" ="precipProbability", "HT" = "temperatureHigh",
    "LT" = "temperatureLow", "DP" = "dewPoint", "HUM" = "humidity", "WS" = "windSpeed", "UVI" = "uvIndex")

geo_root <- geo_root %>%
  rename("LAT" = "GPSlatitude", "LONG" = "GPSlongitude")

prop_root <- prop_root %>%
  rename(
    "LBC" = "Lime_Buffer_Capacity", "pH" = "pH", "WpH" = "Water_pH", "BS" = "Base_Saturation_perc",
    "CEC" = "Cation_Exchange_Capacity", "Ca" = "Ca_ppm", "Cu" = "Cu_ppm", "Fe" = "Fe_ppm", "K" = "K_ppm",
    "Mg" = "Mg_ppm", "Mn" = "Mn_ppm", "Na" = "Na_ppm", "Ni" = "Ni_ppm", "P" = "P_ppm","Pb" = "Pb_ppm", "Zn" = "Zn_ppm", 
    "OM" = "Organic_Matter_perc", "N" = "N_perc", "Mo" = "Mo_ppm", "Cd" = "Cd_ppm", "Cr" = "Cr_ppm")


# Get ASV table
asv_root <- as.data.frame(otu_table(root))
asv_root_t <- t(asv_root)

#Hellinger transformation on species abundance data 
#This turns absolute abundance into relative abundance
#We will call this spec.h
asv_root_t <- decostand(asv_root_t, method = "hellinger")

#check whether to use CCA or RDA
decorana(asv_root_t)

#the length of the first axis (DCA1) is 13.2294, which is greater than 3. 
#This suggests that CCA would be more appropriate


# Check for multicollinearity in environmental variables
# Create a data frame with all environmental variables
env_model_root <- lm(`Growing Degree Days` ~ ., data = env_root)
vif_env_root <- vif(env_model_root)

# Identify variables with VIF > 10
high_vif_env_root <- names(vif_env_root[vif_env_root > 10])

# Check for multicollinearity in root properties
# Due to many variables, use stepwise regression to avoid perfect multicollinearity
prop_model_root <- lm(pH ~ ., data = prop_root)
vif_prop_root <- vif(prop_model_root)
vif_prop_root

# Identify variables with VIF > 10
high_vif_prop_root <- names(vif_prop_root[vif_prop_root > 10])

#Let's get rid of WpH, LT, DP, pH, Ca, , CEC, Mg, PP, OMP 
# Remove high-VIF variables from root properties data frame
prop_root_vif_filtered <- prop_root[, 
                                      !colnames(prop_root) %in% high_vif_prop_root, 
                                      drop = FALSE]

# Check remaining variables
colnames(prop_root_vif_filtered)

# Method 1: CCA analysis using only environmental variables after removing collinear variables
spec.cca_env_root <- cca(asv_root_t ~ ., data = env_root)
print(summary(spec.cca_env_root))

# Method 2: CCA analysis using only root properties
spec.cca_prop_root <- cca(asv_root_t ~ ., data = prop_root_vif_filtered)
print(summary(spec.cca_prop_root))

# Method 3: CCA analysis using only geographic information
spec.cca_geo_root <- cca(asv_root_t ~ ., data = geo_root)
print(summary(spec.cca_geo_root))

# Method 4: Variance partitioning analysis to compare contributions of three variable groups
vp_root <- varpart(asv_root_t, env_root, prop_root, geo_root)
vp_root

# Plot variance partitioning diagram
plot(vp_root, Xnames = c("Environment", "root Properties", "Geography"), bg = c("red", "blue", "green"))


# Use forward selection to select important variables from environmental variables
step.env_root <- ordistep(cca(asv_root_t ~ 1, data = env_root), 
                           scope = formula(spec.cca_env_root), 
                           direction = "forward", 
                           permutations = 999)

# Use forward selection to select important variables from root properties
step.prop_root <- ordistep(cca(asv_root_t ~ 1, data = prop_root), 
                            scope = formula(spec.cca_prop_root), 
                            direction = "forward", 
                            permutations = 999)

# Use forward selection to select important variables from geographic information
step.geo_root <- ordistep(cca(asv_root_t ~ 1, data = geo_root), 
                           scope = formula(spec.cca_geo_root), 
                           direction = "forward", 
                           permutations = 999)

# Extract selected variables
selected_env_vars_root <- attr(terms(step.env_root$call$formula), "term.labels")
selected_prop_vars_root <- attr(terms(step.prop_root$call$formula), "term.labels")
selected_geo_vars_root <- attr(terms(step.geo_root$call$formula), "term.labels")

# Create a data frame with all selected variables
selected_vars_root <- data.frame(
  env_root[, intersect(selected_env_vars_root, colnames(env_root)), drop = FALSE],
  prop_root[, intersect(selected_prop_vars_root, colnames(prop_root)), drop = FALSE],
  geo_root[, intersect(selected_geo_vars_root, colnames(geo_root)), drop = FALSE]
)

# Perform integrated CCA analysis with all selected variables
spec.cca_combined_root <- cca(asv_root_t ~ ., data = selected_vars_root)
print(summary(spec.cca_combined_root))

# Extract CCA results
cca_combined_root <- scores(spec.cca_combined_root, display = c("sites", "species", "bp"), scaling = 1)

# Extract site and environmental variable coordinates
sites_combined_root <- as.data.frame(cca_combined_root$sites)
bp_combined_root <- as.data.frame(cca_combined_root$biplot)

# Merge site coordinates with metadata
samples_combined_root <- cbind(sites_combined_root, sample_data_root_df)

# Calculate percentage of variance explained by CCA axes
aa_combined_root <- summary(spec.cca_combined_root)$cont$importance[2, 1:2]
aa_combined_root <- rbind(aa_combined_root, aa_combined_root * 100)

# Create environmental variables data frame for plotting
envi_combined_root <- data.frame(
  env = rownames(bp_combined_root),
  CCA1 = bp_combined_root$CCA1,
  CCA2 = bp_combined_root$CCA2
)

# Add variable type colors
envi_combined_root$type <- "Other"
envi_combined_root$type[envi_combined_root$env %in% selected_env_vars_root] <- "Environment"
envi_combined_root$type[envi_combined_root$env %in% selected_prop_vars_root] <- "root"
envi_combined_root$type[envi_combined_root$env %in% selected_geo_vars_root] <- "Geography"

# Plot integrated CCA diagram
cca_root <- ggplot() + 
  geom_segment(data = envi_combined_root, aes(x = 0, y = 0, xend = CCA1*3.1, yend = CCA2*3.1, color = type), 
               arrow = arrow(length = unit(0.5, "cm")), linewidth = 1.2) + 
  geom_text(data = envi_combined_root, aes(x = CCA1*3.1, y = CCA2*3.1, label = env, color = type), 
            size = 10, check_overlap = FALSE) +
  geom_point(data = samples_combined_root, aes(x=CCA1, y=CCA2, shape = factor(`Growing_Degree_Days`)), size = 7) +
  scale_color_manual(values = c("Environment" = "red", "root" = "blue", "Geography" = "green", "Other" = "gray")) +
  theme_bw() + 
  theme(panel.grid.major=element_line(colour=NA), panel.grid.minor = element_blank()) + 
  xlab(paste('CCA1 (', aa_combined_root[2,1], '%)', sep ='')) + 
  ylab(paste('CCA2 (', aa_combined_root[2,2], '%)', sep ='')) + 
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", face = "bold", size = 15),      
        axis.text.x=element_text(angle=0,hjust=1,vjust=1.0,size=15, colour = "black"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black"),
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 10, face = "bold", angle = 0),
        legend.position= "right",
        legend.title = element_text(face="bold", size = 15)) +
  labs(color = "Variable Type", shape = "Growing Degree Days") +
  ggtitle("root Microbiome CCA Analysis")

cca_root

#save plot in pdf and png
ggsave("cca_root.png", plot = cca_root , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("cca_root.pdf", plot = cca_root, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")


# Test integrated model significance
# ANOVA - test significance of the entire CCA model
anova_combined_root <- anova(spec.cca_combined_root, perm.max=999)

# Test significance of each variable
anova_terms_root <- anova(spec.cca_combined_root, by="terms", perm.max=999)

# Extract raw p-values excluding NA values
raw_pvals_root <- anova_terms_root$`Pr(>F)`[!is.na(anova_terms_root$`Pr(>F)`)]

# Apply Bonferroni correction
bonf_pvals_root <- p.adjust(raw_pvals_root, method = "bonferroni")

# Apply other correction methods for comparison
fdr_pvals_root <- p.adjust(raw_pvals_root, method = "fdr")
holm_pvals_root <- p.adjust(raw_pvals_root, method = "holm")

# Create complete correction results table
correction_results_root <- data.frame(
  Variable = rownames(anova_terms_root)[!is.na(anova_terms_root$`Pr(>F)`)],
  F_value = anova_terms_root$F[!is.na(anova_terms_root$`Pr(>F)`)],
  Raw_P_value = raw_pvals_root,
  Bonferroni_P_value = bonf_pvals_root,
  FDR_P_value = fdr_pvals_root,
  Holm_P_value = holm_pvals_root
)

# Add variable type information
correction_results_root$Type <- "Other"
correction_results_root$Type[correction_results_root$Variable %in% selected_env_vars_root] <- "Environment"
correction_results_root$Type[correction_results_root$Variable %in% selected_prop_vars_root] <- "root"
correction_results_root$Type[correction_results_root$Variable %in% selected_geo_vars_root] <- "Geography"

# Calculate adjusted R-squared
rsq_root <- RsquareAdj(spec.cca_combined_root)

# Save results table to CSV file
write.csv(correction_results_root, "C:/Aduragbemi/Manuscript/Review 5/PDF5/cca_correction_results_root.csv", row.names = FALSE)




