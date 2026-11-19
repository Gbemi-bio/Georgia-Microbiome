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


#=============================================
# Alpha diversity vs environmental variables 
# Multiple-testing correction using BH-FDR
#=============================================

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


#===========================================
# Define predictor groups
#===========================================

climate_vars <- c(
  "temperatureHigh",
  "temperatureLow",
  "precipIntensity",
  "windSpeed",
  "humidity",
  "uvIndex",
  "dewPoint"
)

geography_vars <- c(
  "GPSlatitude",
  "GPSlongitude"
)

soil_vars <- c(
  "Lime_Buffer_Capacity",
  "pH",
  "Water_pH",
  "Base_Saturation_perc",
  "Cation_Exchange_Capacity",
  "Ca_ppm",
  "K_ppm",
  "Mg_ppm",
  "Na_ppm",
  "P_ppm",
  "Cd_ppm",
  "Cr_ppm",
  "Cu_ppm",
  "Fe_ppm",
  "Mn_ppm",
  "Mo_ppm",
  "Ni_ppm",
  "Pb_ppm",
  "Zn_ppm",
  "Organic_Matter_perc",
  "N_perc"
)

# Alpha diversity metrics
metrics <- c("Shannon","PD")

#===========================================
# Function to run GLMs
#===========================================

run_glm_group <- function(data, predictors, group_name){
  
  results <- data.frame()
  
  for(metric in metrics){
    
    for(var in predictors){
      
      model <- glm(
        as.formula(paste(metric,"~",var)),
        data = data
      )
      
      coef_table <- summary(model)$coefficients
      
      results <- rbind(
        results,
        data.frame(
          Group = group_name,
          Response = metric,
          Predictor = var,
          Estimate = coef_table[2,1],
          Std_Error = coef_table[2,2],
          t_value = coef_table[2,3],
          P_value = coef_table[2,4]
        )
      )
      
    }
    
  }
  
  ## BH correction WITHIN each response variable
  results$FDR <- NA
  
  for(m in unique(results$Response)){
    
    idx <- results$Response == m
    
    results$FDR[idx] <-
      p.adjust(results$P_value[idx],
               method="BH")
    
  }
  
  results
}

#===========================
# Climate 
#===========================

results_climate <-
  run_glm_group(
    soil_continuous_alpha_1,
    climate_vars,
    "Climate"
  )

results_climate

#=================================
# Geography
#=================================

results_geography <-
  run_glm_group(
    soil_continuous_alpha_1,
    geography_vars,
    "Geography"
  )

results_geography

#===========================================
# Soil properties
#===========================================

results_soilchem <-
  run_glm_group(
    soil_continuous_alpha_1,
    soil_vars,
    "Soil Properties"
  )

results_soilchem

#save each table
write.csv(
  results_climate,
  "C:/Aduragbemi/Manuscript/Review 7/PDF7/Soil_Climate_GLM_FDR.csv",
  row.names=FALSE
)

write.csv(
  results_geography,
  "C:/Aduragbemi/Manuscript/Review 7/PDF7/Soil_Geography_GLM_FDR.csv",
  row.names=FALSE
)

write.csv(
  results_soilchem,
  "C:/Aduragbemi/Manuscript/Review 7/PDF7/Soil_SoilProperties_GLM_FDR.csv",
  row.names=FALSE
)


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


#===========================================
# Define predictor groups
#===========================================

climate_vars <- c(
  "temperatureHigh",
  "temperatureLow",
  "precipIntensity",
  "windSpeed",
  "humidity",
  "uvIndex",
  "dewPoint"
)

geography_vars <- c(
  "GPSlatitude",
  "GPSlongitude"
)

rhizo_vars <- c(
  "Lime_Buffer_Capacity",
  "pH",
  "Water_pH",
  "Base_Saturation_perc",
  "Cation_Exchange_Capacity",
  "Ca_ppm",
  "K_ppm",
  "Mg_ppm",
  "Na_ppm",
  "P_ppm",
  "Cd_ppm",
  "Cr_ppm",
  "Cu_ppm",
  "Fe_ppm",
  "Mn_ppm",
  "Mo_ppm",
  "Ni_ppm",
  "Pb_ppm",
  "Zn_ppm",
  "Organic_Matter_perc",
  "N_perc"
)

# Alpha diversity metrics
metrics <- c("Shannon","PD")

#===========================================
# Function to run GLMs
#===========================================

run_glm_group <- function(data, predictors, group_name){
  
  results <- data.frame()
  
  for(metric in metrics){
    
    for(var in predictors){
      
      model <- glm(
        as.formula(paste(metric,"~",var)),
        data = data
      )
      
      coef_table <- summary(model)$coefficients
      
      results <- rbind(
        results,
        data.frame(
          Group = group_name,
          Response = metric,
          Predictor = var,
          Estimate = coef_table[2,1],
          Std_Error = coef_table[2,2],
          t_value = coef_table[2,3],
          P_value = coef_table[2,4]
        )
      )
      
    }
    
  }
  
  ## BH correction WITHIN each response variable
  results$FDR <- NA
  
  for(m in unique(results$Response)){
    
    idx <- results$Response == m
    
    results$FDR[idx] <-
      p.adjust(results$P_value[idx],
               method="BH")
    
  }
  
  results
}

#===========================
# Climate 
#===========================

results_climate <-
  run_glm_group(
    rhizo_continuous_alpha_1,
    climate_vars,
    "Climate"
  )

results_climate

#=================================
# Geography
#=================================

results_geography <-
  run_glm_group(
    rhizo_continuous_alpha_1,
    geography_vars,
    "Geography"
  )

results_geography

#===========================================
# rhizo properties
#===========================================

results_rhizochem <-
  run_glm_group(
    rhizo_continuous_alpha_1,
    rhizo_vars,
    "rhizo Properties"
  )

results_rhizochem

#save each table
write.csv(
  results_climate,
  "C:/Aduragbemi/Manuscript/Review 7/PDF7/rhizo_Climate_GLM_FDR.csv",
  row.names=FALSE
)

write.csv(
  results_geography,
  "C:/Aduragbemi/Manuscript/Review 7/PDF7/rhizo_Geography_GLM_FDR.csv",
  row.names=FALSE
)

write.csv(
  results_rhizochem,
  "C:/Aduragbemi/Manuscript/Review 7/PDF7/rhizo_rhizoProperties_GLM_FDR.csv",
  row.names=FALSE
)



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


#===========================================
# Define predictor groups
#===========================================

climate_vars <- c(
  "temperatureHigh",
  "temperatureLow",
  "precipIntensity",
  "windSpeed",
  "humidity",
  "uvIndex",
  "dewPoint"
)

geography_vars <- c(
  "GPSlatitude",
  "GPSlongitude"
)

root_vars <- c(
  "Lime_Buffer_Capacity",
  "pH",
  "Water_pH",
  "Base_Saturation_perc",
  "Cation_Exchange_Capacity",
  "Ca_ppm",
  "K_ppm",
  "Mg_ppm",
  "Na_ppm",
  "P_ppm",
  "Cd_ppm",
  "Cr_ppm",
  "Cu_ppm",
  "Fe_ppm",
  "Mn_ppm",
  "Mo_ppm",
  "Ni_ppm",
  "Pb_ppm",
  "Zn_ppm",
  "Organic_Matter_perc",
  "N_perc"
)

# Alpha diversity metrics
metrics <- c("Shannon","PD")

#===========================================
# Function to run GLMs
#===========================================

run_glm_group <- function(data, predictors, group_name){
  
  results <- data.frame()
  
  for(metric in metrics){
    
    for(var in predictors){
      
      model <- glm(
        as.formula(paste(metric,"~",var)),
        data = data
      )
      
      coef_table <- summary(model)$coefficients
      
      results <- rbind(
        results,
        data.frame(
          Group = group_name,
          Response = metric,
          Predictor = var,
          Estimate = coef_table[2,1],
          Std_Error = coef_table[2,2],
          t_value = coef_table[2,3],
          P_value = coef_table[2,4]
        )
      )
      
    }
    
  }
  
  ## BH correction WITHIN each response variable
  results$FDR <- NA
  
  for(m in unique(results$Response)){
    
    idx <- results$Response == m
    
    results$FDR[idx] <-
      p.adjust(results$P_value[idx],
               method="BH")
    
  }
  
  results
}

#===========================
# Climate 
#===========================

results_climate <-
  run_glm_group(
    root_continuous_alpha_1,
    climate_vars,
    "Climate"
  )

results_climate

#=================================
# Geography
#=================================

results_geography <-
  run_glm_group(
    root_continuous_alpha_1,
    geography_vars,
    "Geography"
  )

results_geography

#===========================================
# root properties
#===========================================

results_rootchem <-
  run_glm_group(
    root_continuous_alpha_1,
    root_vars,
    "root Properties"
  )

results_rootchem

#save each table
write.csv(
  results_climate,
  "C:/Aduragbemi/Manuscript/Review 7/PDF7/root_Climate_GLM_FDR.csv",
  row.names=FALSE
)

write.csv(
  results_geography,
  "C:/Aduragbemi/Manuscript/Review 7/PDF7/root_Geography_GLM_FDR.csv",
  row.names=FALSE
)

write.csv(
  results_rootchem,
  "C:/Aduragbemi/Manuscript/Review 7/PDF7/root_rootProperties_GLM_FDR.csv",
  row.names=FALSE
)



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


#===========================================
# Define predictor groups
#===========================================

climate_vars <- c(
  "temperatureHigh",
  "temperatureLow",
  "precipIntensity",
  "windSpeed",
  "humidity",
  "uvIndex",
  "dewPoint"
)

geography_vars <- c(
  "GPSlatitude",
  "GPSlongitude"
)

stalk_vars <- c(
  "Lime_Buffer_Capacity",
  "pH",
  "Water_pH",
  "Base_Saturation_perc",
  "Cation_Exchange_Capacity",
  "Ca_ppm",
  "K_ppm",
  "Mg_ppm",
  "Na_ppm",
  "P_ppm",
  "Cd_ppm",
  "Cr_ppm",
  "Cu_ppm",
  "Fe_ppm",
  "Mn_ppm",
  "Mo_ppm",
  "Ni_ppm",
  "Pb_ppm",
  "Zn_ppm",
  "Organic_Matter_perc",
  "N_perc"
)

# Alpha diversity metrics
metrics <- c("Shannon","PD")

#===========================================
# Function to run GLMs
#===========================================

run_glm_group <- function(data, predictors, group_name){
  
  results <- data.frame()
  
  for(metric in metrics){
    
    for(var in predictors){
      
      model <- glm(
        as.formula(paste(metric,"~",var)),
        data = data
      )
      
      coef_table <- summary(model)$coefficients
      
      results <- rbind(
        results,
        data.frame(
          Group = group_name,
          Response = metric,
          Predictor = var,
          Estimate = coef_table[2,1],
          Std_Error = coef_table[2,2],
          t_value = coef_table[2,3],
          P_value = coef_table[2,4]
        )
      )
      
    }
    
  }
  
  ## BH correction WITHIN each response variable
  results$FDR <- NA
  
  for(m in unique(results$Response)){
    
    idx <- results$Response == m
    
    results$FDR[idx] <-
      p.adjust(results$P_value[idx],
               method="BH")
    
  }
  
  results
}

#===========================
# Climate 
#===========================

results_climate <-
  run_glm_group(
    stalk_continuous_alpha_1,
    climate_vars,
    "Climate"
  )

results_climate

#=================================
# Geography
#=================================

results_geography <-
  run_glm_group(
    stalk_continuous_alpha_1,
    geography_vars,
    "Geography"
  )

results_geography

#===========================================
# stalk properties
#===========================================

results_stalkchem <-
  run_glm_group(
    stalk_continuous_alpha_1,
    stalk_vars,
    "stalk Properties"
  )

results_stalkchem

#save each table
write.csv(
  results_climate,
  "C:/Aduragbemi/Manuscript/Review 7/PDF7/stalk_Climate_GLM_FDR.csv",
  row.names=FALSE
)

write.csv(
  results_geography,
  "C:/Aduragbemi/Manuscript/Review 7/PDF7/stalk_Geography_GLM_FDR.csv",
  row.names=FALSE
)

write.csv(
  results_stalkchem,
  "C:/Aduragbemi/Manuscript/Review 7/PDF7/stalk_stalkProperties_GLM_FDR.csv",
  row.names=FALSE
)