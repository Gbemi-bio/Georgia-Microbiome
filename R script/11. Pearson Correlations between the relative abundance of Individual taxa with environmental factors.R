# Load required libraries
library(phyloseq)
library(microeco)
library(ggplot2)
library(file2meco)

set.seed(1234)

# =======================================================================
# Correlation between environmental factors and individual class groups
# =======================================================================

#for soil samples
soil_bac_sub_no_bad_Filtered_final_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/soil_bac_sub_no_bad_Filtered_final_non_normalized.rds")

soil_norm = transform_sample_counts(soil_bac_sub_no_bad_Filtered_final_non_normalized, function(x) x / sum(x) )
soil_network = filter_taxa(soil_norm, function(x) mean(x) > 0.0001, TRUE)
soil_network

#first convert to microtable
soil_amp <- phyloseq2meco(soil_network)
soil_amp

soil_env <- soil_amp[["sample_table"]][, 53:60]

t1_soil <- trans_env$new(dataset = soil_amp, add_data = soil_env)

# 'p_adjust_type = "Env"' means p adjustment is performed for each environmental variable separately.
t1_soil$cal_cor(use_data = "Class", by_group = "Growing_Degree_Days", p_adjust_method = "fdr", p_adjust_type = "Env", cor_method = "pearson")


t1_soil$res_cor

# default ggplot2 method with clustering
soil_cor <- t1_soil$plot_cor()
soil_cor

#change the order of some variables
soil_cor$data$by_group <- factor(soil_cor$data$by_group, levels=c("600", "1400"))
soil_cor$data$by_group
soil_cor

# include the class that show only significant with climatic variables
soil_cor$data <- soil_cor$data[
  soil_cor$data$Taxa %in% c(
    "c__Vicinamibacteria",
    "c__Verrucomicrobiae",
    "c__TK10",
    "c__Thermoleophilia",
    "c__Planctomycetes",
    "c__Phycisphaerae",
    "c__Nitrospiria",
    "c__Negativicutes",
    "c__Methylomirabilia",
    "c__MB-A2-108",
    "c__Latescibacterota",
    "c__Ktedonobacteria",
    "c__KD4-96",
    "c__Holophagae",
    "c__Gemmatimonadetes",
    "c__Gammaproteobacteria",
    "c__Chloroflexia",
    "c__Blastocatellia",
    "c__Berkelbacteria",
    "c__Bacilli",
    "c__Anaerolineae",
    "c__Anaproteobacteria",
    "c__Actinobacteria",
    "c__Acidobacteriae",
    "c__Acidimicrobiia"
  ),
]

#rename the taxaby removing the C__
soil_cor$data$Taxa <- gsub("^c__", "", soil_cor$data$Taxa)
soil_cor

#include only the environmental variables that show significant correlations with class
soil_cor$data <- soil_cor$data[(soil_cor$data$Env %in% c(
  "precipIntensity", 
  "precipProbability"
)), 
]
soil_cor


soil_cor$data$Env <- factor(soil_cor$data$Env, labels = c(
  "precipIntensity" = "Precipitation Intensity",
  "precipProbability" = "Precipitation Probability"
))
soil_cor

soil_cor_final <- soil_cor + theme(plot.title = element_text(family = "serif", face = "bold", size = (70)),
                                   axis.text.x=element_text(angle=45,hjust=1,vjust=1.0,size=60, colour = "black"),
                                   axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=80, colour = "black"),
                                   legend.text = element_text(size = 45, face = "bold"),
                                   strip.text = element_text(size = 70, face = "bold", angle = 0),
                                   legend.title = element_text(face="bold", size = 70),
                                   legend.key.height = unit(5, 'cm'), legend.key.width = unit(5, 'cm')) +
                                   ggtitle("soil")


soil_cor_final

#save plot
ggsave("soil_cor_final.png", plot =soil_cor_final, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 30, height = 30, units = c("in"), device = "png")

ggsave("soil_cor_final.pdf", plot = soil_cor_final, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 30, height = 30, units = c("in"), device = "pdf")




#for Rhizosphere samples
rhizo_bac_sub_no_bad_Filtered_final_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/rhizo_bac_sub_no_bad_Filtered_final_non_normalized.rds")

rhizo_norm = transform_sample_counts(rhizo_bac_sub_no_bad_Filtered_final_non_normalized, function(x) x / sum(x) )
rhizo_network = filter_taxa(rhizo_norm, function(x) mean(x) > 0.0001, TRUE)
rhizo_network

#first convert to microtable
rhizo_amp <- phyloseq2meco(rhizo_network)
rhizo_amp

rhizo_env <- rhizo_amp[["sample_table"]][, 53:60]

t1_rhizo <- trans_env$new(dataset = rhizo_amp, add_data = rhizo_env)

# 'p_adjust_type = "Env"' means p adjustment is performed for each environmental variable separately.
t1_rhizo$cal_cor(use_data = "Class", by_group = "Growing_Degree_Days", p_adjust_method = "fdr", p_adjust_type = "Env", cor_method = "pearson")


t1_rhizo$res_cor

# default ggplot2 method with clustering
rhizo_cor <- t1_rhizo$plot_cor()
rhizo_cor

#change the order of some variables
rhizo_cor$data$by_group <- factor(rhizo_cor$data$by_group, levels=c("600", "1400"))
rhizo_cor$data$by_group
rhizo_cor

# include the class that show only significant with climatic variables
rhizo_cor$data <- rhizo_cor$data[
  rhizo_cor$data$Taxa %in% c(
    "c__Thermoleophilia",
    "c__Phycisphaerae",
    "c__Gemmatimonadetes",
    "c__Gammaproteobacteria",
    "c__Chloroflexia",
    "c__Blastocatellia",
    "c__Bacilli",
    "c__Alphaproteobacteria",
    "c__Actinobacteria",
    "c__Acidobacteriae"
  ),
]
rhizo_cor

#rename the taxaby removing the C__
rhizo_cor$data$Taxa <- gsub("^c__", "", rhizo_cor$data$Taxa)
rhizo_cor


#include only the environmental variables that show significant correlations with class
rhizo_cor$data <- rhizo_cor$data[(rhizo_cor$data$Env %in% c(
  "precipIntensity", 
  "precipProbability", 
  "uvIndex"
)), 
]
rhizo_cor


rhizo_cor$data$Env <- factor(rhizo_cor$data$Env, labels = c(
  "precipIntensity" = "Precipitation Intensity",
  "precipProbability" = "Precipitation Probability",
  "uvIndex" = "Ultraviolet Index"
))
rhizo_cor

rhizo_cor_final <- rhizo_cor + theme(plot.title = element_text(family = "serif", face = "bold", size = (70)),
                                     axis.text.x=element_text(angle=45,hjust=1,vjust=1.0,size=60, colour = "black"),
                                     axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=80, colour = "black"),
                                     legend.text = element_text(size = 45, face = "bold"),
                                     strip.text = element_text(size = 70, face = "bold", angle = 0),
                                     legend.title = element_text(face="bold", size = 70),
                                     legend.key.height = unit(5, 'cm'), legend.key.width = unit(5, 'cm')) +
                                     ggtitle("Rhizosphere")


rhizo_cor_final

#save plot
ggsave("rhizo_cor_final.png", plot =rhizo_cor_final, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 30, height = 30, units = c("in"), device = "png")

ggsave("rhizo_cor_final.pdf", plot = rhizo_cor_final, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 30, height = 30, units = c("in"), device = "pdf")



#for Root samples
root_bac_sub_no_bad_Filtered_final_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/root_bac_sub_no_bad_Filtered_final_non_normalized.rds")

root_norm = transform_sample_counts(root_bac_sub_no_bad_Filtered_final_non_normalized, function(x) x / sum(x) )
root_network = filter_taxa(root_norm, function(x) mean(x) > 0.0001, TRUE)
root_network

#first convert to microtable
root_amp <- phyloseq2meco(root_network)
root_amp

root_env <- root_amp[["sample_table"]][, 53:60]
root_env

t1_root <- trans_env$new(dataset = root_amp, add_data = root_env)

# 'p_adjust_type = "Env"' means p adjustment is performed for each environmental variable separately.
t1_root$cal_cor(use_data = "Class", by_group = "Growing_Degree_Days", p_adjust_method = "fdr", p_adjust_type = "Env", cor_method = "pearson")


t1_root$res_cor

# default ggplot2 method with clustering
root_cor <- t1_root$plot_cor()
root_cor

#change the order of some variables
root_cor$data$by_group <- factor(root_cor$data$by_group, levels=c("600", "1400"))
root_cor$data$by_group
root_cor

# include the class that show only significant with climatic variables
root_cor$data <- root_cor$data[
  root_cor$data$Taxa %in% c(
    "c__Chloroflexia",
    "c__Bacteroidia",
    "c__Actinobacteria"
  ),
]
root_cor

#rename the taxaby removing the C__
root_cor$data$Taxa <- gsub("^c__", "", root_cor$data$Taxa)
root_cor

#include only the environmental variables that show significant correlations with class
root_cor$data <- root_cor$data[(root_cor$data$Env %in% c("precipIntensity", "precipProbability")), ]
root_cor


root_cor$data$Env <- factor(root_cor$data$Env, labels = c("precipIntensity" = "Precipitation Intensity", "precipProbability"= "Precipitation Probability"))
root_cor

root_cor_final <- root_cor + theme(plot.title = element_text(family = "serif", face = "bold", size = (70)),
                                   axis.text.x=element_text(angle=45,hjust=1,vjust=1.0,size=60, colour = "black"),
                                   axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=80, colour = "black"),
                                   legend.text = element_text(size = 45, face = "bold"),
                                   strip.text = element_text(size = 70, face = "bold", angle = 0),
                                   legend.title = element_text(face="bold", size = 70),
                                   legend.key.height = unit(5, 'cm'), legend.key.width = unit(5, 'cm')) +
                                   ggtitle("Root")


root_cor_final

#save plot
ggsave("root_cor_final.png", plot =root_cor_final, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 30, height = 30, units = c("in"), device = "png")

ggsave("root_cor_final.pdf", plot = root_cor_final, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 30, height = 30, units = c("in"), device = "pdf")






#for stalk samples
stalk_bac_sub_no_bad_Filtered_final_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/stalk_bac_sub_no_bad_Filtered_final_non_normalized.rds")

stalk_norm = transform_sample_counts(stalk_bac_sub_no_bad_Filtered_final_non_normalized, function(x) x / sum(x) )
stalk_network = filter_taxa(stalk_norm, function(x) mean(x) > 0.0001, TRUE)
stalk_network

#first convert to microtable
stalk_amp <- phyloseq2meco(stalk_network)
stalk_amp

stalk_env <- stalk_amp[["sample_table"]][, 53:60]

t1_stalk <- trans_env$new(dataset = stalk_amp, add_data = stalk_env)

# 'p_adjust_type = "Env"' means p adjustment is performed for each environmental variable separately.
t1_stalk$cal_cor(use_data = "Class", by_group = "Growing_Degree_Days", p_adjust_method = "fdr", p_adjust_type = "Env", cor_method = "pearson")
t1_stalk$res_cor

#default ggplot2 method with clustering
stalk_cor <- t1_stalk$plot_cor()
stalk_cor

#change the order of some variables
stalk_cor$data$by_group <- factor(stalk_cor$data$by_group, levels=c("600", "1400"))
stalk_cor$data$by_group
stalk_cor


# include the class that show only significant with climatic variables
stalk_cor$data <- stalk_cor$data[
  (stalk_cor$data$Taxa %in% c(
    "c__Planctomycetes",
    "c__Gammaproteobacteria",
    "c__Chloroflexia",
    "c__Bdellovibrionia",
    "c__Bacteroidia",
    "c__Armatimonadia",
    "c__Alphaproteobacteria",
    "c__Actinobacteria"
  )),
]

#rename the taxaby removing the C__
stalk_cor$data$Taxa <- gsub("^c__", "", stalk_cor$data$Taxa)
stalk_cor

#include only the environmental variables that show significant correlations with class
stalk_cor$data <- stalk_cor$data[(stalk_cor$data$Env %in% c(
  "precipIntensity", 
  "precipProbability"
)), 
]
stalk_cor


stalk_cor$data$Env <- factor(stalk_cor$data$Env, labels = c(
  "precipIntensity" = "Precipitation Intensity",
  "precipProbability" = "Precipitation Probability"
))
stalk_cor

stalk_cor_final <- stalk_cor + theme(plot.title = element_text(family = "serif", face = "bold", size = (70)),
                                     axis.text.x=element_text(angle=45,hjust=1,vjust=1.0,size=60, colour = "black"),
                                     axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=80, colour = "black"),
                                     legend.text = element_text(size = 45, face = "bold"),
                                     strip.text = element_text(size = 70, face = "bold", angle = 0),
                                     legend.title = element_text(face="bold", size = 70),
                                     legend.key.height = unit(5, 'cm'), legend.key.width = unit(5, 'cm')) +
                                     ggtitle("stalk")


stalk_cor_final

#save plot
ggsave("stalk_cor_final.png", plot =stalk_cor_final, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 30, height = 30, units = c("in"), device = "png")

ggsave("stalk_cor_final.pdf", plot = stalk_cor_final, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 30, height = 30, units = c("in"), device = "pdf")





# =======================================================================
# Correlation between soil chemistry and individual class groups
# =======================================================================

#For soil samples

soil_bac_sub_no_bad_Filtered_final_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/soil_bac_sub_no_bad_Filtered_final_non_normalized.rds")

soil_norm = transform_sample_counts(soil_bac_sub_no_bad_Filtered_final_non_normalized, function(x) x / sum(x) )
soil_network = filter_taxa(soil_norm, function(x) mean(x) > 0.0001, TRUE)
soil_network

#first convert to microtable
soil_amp <- phyloseq2meco(soil_network)
soil_amp

soil_env <- soil_amp[["sample_table"]][, 26:46]
soil_env
t1_soil <- trans_env$new(dataset = soil_amp, add_data = soil_env)

# 'p_adjust_type = "Env"' means p adjustment is performed for each environmental variable separately.
t1_soil$cal_cor(use_data = "Class", by_group = "Growing_Degree_Days", p_adjust_method = "fdr", p_adjust_type = "Env", cor_method = "pearson")

t1_soil$res_cor

# default ggplot2 method with clustering
soil_prop <- t1_soil$plot_cor()
soil_prop

soil_prop[["data"]][["by_group"]] <- factor(soil_prop[["data"]][["by_group"]], levels=c("600", "1400"))
soil_prop



soil_prop$data <- soil_prop$data[
  (soil_prop$data$Taxa %in% c(
    "c__Vicinamibacteria",
    "c__Verrucomicrobiae",
    "c__TK10",
    "c__Thermoleophilia",
    "c__Planctomycetes",
    "c__Phycisphaerae",
    "c__Nitrospiria",
    "c__Myxococcia",
    "c__Methylomirabilia",
    "c__MB-A2-108",
    "c__Latescibacterota",
    "c__Ktedonobacteria",
    "c__KD4-96",
    "c__JG30-KF-CM66",
    "c__Holophagae",
    "c__Gitt-GS-136",
    "c__Gemmatimonadetes",
    "c__Gammaproteobacteria",
    "c__Chloroflexia",
    "c__Blastocatellia",
    "c__Bacteroidia",
    "c__bacteria25",
    "c__Bacilli",
    "c__Babeliae",
    "c__Anaerolineae",
    "c__Alphaproteobacteria",
    "c__AD3",
    "c__Actinobacteria",
    "c__Acidobacteriae",
    "c__Acidimicrobiia"
  )),
]

soil_prop

#rename the taxaby removing the C__
soil_prop$data$Taxa <- gsub("^c__", "", soil_prop$data$Taxa)
soil_prop

#include only the soil chemistry variables that show significant correlations with class
soil_prop$data <- soil_prop$data[(soil_prop$data$Env %in% c(
  "Ca_ppm", 
  "Cation_Exchange_Capacity",
  "Cr_ppm", 
  "Fe_ppm", 
  "K_ppm",
  "Lime_Buffer_Capacity",
  "Mn_ppm",
  "Na_ppm", 
  "P_ppm"
)), 
]
soil_prop

soil_prop$data$Env <- factor(soil_prop$data$Env, labels = c(
  "Ca_ppm" = "Calcium",
  "Cation_Exchange_Capacity" = "Cation Exchange Capacity",
  "Cr_ppm" = "Chromium",
  "Fe_ppm" = "Iron",
  "K_ppm" = "Potassium",
  "Lime_Buffer_Capacity" = "Lime Buffer Capacity",
  "Mn_ppm" = "Manganese",
  "Na_ppm" = "Sodium",
  "P_ppm" = "Phosphorus"
))
soil_prop


soil_prop_final <- soil_prop + theme(plot.title = element_text(family = "serif", face = "bold", size = (80)),
                                     axis.text.x=element_text(angle=45,hjust=1,vjust=1.0,size=50, colour = "black"),
                                     axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=80, colour = "black"),
                                     legend.text = element_text(size = 45, face = "bold"),
                                     strip.text = element_text(size = 70, face = "bold", angle = 0),
                                     legend.title = element_text(face="bold", size = 70),
                                     legend.key.height = unit(5, 'cm'), legend.key.width = unit(5, 'cm')) +
                                     ggtitle("Soil")


soil_prop_final



#save plot
ggsave("soil_prop_final.png", plot =soil_prop_final, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 30, height = 30, units = c("in"), device = "png")

ggsave("soil_prop_final.pdf", plot = soil_prop_final, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 30, height = 30, units = c("in"), device = "pdf")




#For Rhizosphere samples

rhizo_bac_sub_no_bad_Filtered_final_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/rhizo_bac_sub_no_bad_Filtered_final_non_normalized.rds")

rhizo_norm = transform_sample_counts(rhizo_bac_sub_no_bad_Filtered_final_non_normalized, function(x) x / sum(x) )
rhizo_network = filter_taxa(rhizo_norm, function(x) mean(x) > 0.0001, TRUE)
rhizo_network

#first convert to microtable
rhizo_amp <- phyloseq2meco(rhizo_network)
rhizo_amp

rhizo_env <- rhizo_amp[["sample_table"]][, 26:46]
rhizo_env
t1_rhizo <- trans_env$new(dataset = rhizo_amp, add_data = rhizo_env)

# 'p_adjust_type = "Env"' means p adjustment is performed for each environmental variable separately.
t1_rhizo$cal_cor(use_data = "Class", by_group = "Growing_Degree_Days", p_adjust_method = "fdr", p_adjust_type = "Env", cor_method = "pearson")

t1_rhizo$res_cor

# default ggplot2 method with clustering
rhizo_prop <- t1_rhizo$plot_cor()
rhizo_prop

rhizo_prop[["data"]][["by_group"]] <- factor(rhizo_prop[["data"]][["by_group"]], levels=c("600", "1400"))
rhizo_prop



rhizo_prop$data <- rhizo_prop$data[
  (rhizo_prop$data$Taxa %in% c(
    "c__Vicinamibacteria",
    "c__Verrucomicrobiae",
    "c__TK10",
    "c__Thermoleophilia",
    "c__Planctomycetes",
    "c__Phycisphaerae",
    "c__OLB14",
    "c__Myxococcia",
    "c__Methylomirabilia",
    "c__Ktedonobacteria",
    "c__KD4-96",
    "c__Gitt-GS-136",
    "c__Gemmatimonadetes",
    "c__Gammaproteobacteria",
    "c__Chloroflexia",
    "c__Blastocatellia",
    "c__Bacteroidia",
    "c__Bacilli",
    "c__Anaerolineae",
    "c__Alphaproteobacteria",
    "c__AD3",
    "c__Actinobacteria",
    "c__Acidobacteriae"
  )),
]


rhizo_prop

#rename the taxaby removing the C__
rhizo_prop$data$Taxa <- gsub("^c__", "", rhizo_prop$data$Taxa)
rhizo_prop


#include only the physicochemical variables that show significant correlations with class
rhizo_prop$data <- rhizo_prop$data[(rhizo_prop$data$Env %in% c(
  "Ca_ppm", 
  "Cation_Exchange_Capacity",
  "Cr_ppm", 
  "Fe_ppm", 
  "K_ppm",
  "Lime_Buffer_Capacity",
  "Mn_ppm",
  "Na_ppm", 
  "P_ppm"
)), 
]
rhizo_prop

rhizo_prop$data$Env <- factor(rhizo_prop$data$Env, labels = c(
  "Ca_ppm" = "Calcium",
  "Cation_Exchange_Capacity" = "Cation Exchange Capacity",
  "Cr_ppm" = "Chromium",
  "Fe_ppm" = "Iron",
  "K_ppm" = "Potassium",
  "Lime_Buffer_Capacity" = "Lime Buffer Capacity",
  "Mn_ppm" = "Manganese",
  "Na_ppm" = "Sodium",
  "P_ppm" = "Phosphorus"
))

rhizo_prop


rhizo_prop_final <- rhizo_prop + theme(plot.title = element_text(family = "serif", face = "bold", size = (80)),
                                       axis.text.x=element_text(angle=45,hjust=1,vjust=1.0,size=50, colour = "black"),
                                       axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=80, colour = "black"),
                                       legend.text = element_text(size = 45, face = "bold"),
                                       strip.text = element_text(size = 70, face = "bold", angle = 0),
                                       legend.title = element_text(face="bold", size = 70),
                                       legend.key.height = unit(5, 'cm'), legend.key.width = unit(5, 'cm')) +
                                       ggtitle("Rhizosphere")




rhizo_prop_final


#save plot
ggsave("rhizo_prop_final.png", plot =rhizo_prop_final, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 30, height = 30, units = c("in"), device = "png")

ggsave("rhizo_prop_final.pdf", plot = rhizo_prop_final, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 30, height = 30, units = c("in"), device = "pdf")





#For root samples

root_bac_sub_no_bad_Filtered_final_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/root_bac_sub_no_bad_Filtered_final_non_normalized.rds")

root_norm = transform_sample_counts(root_bac_sub_no_bad_Filtered_final_non_normalized, function(x) x / sum(x) )
root_network = filter_taxa(root_norm, function(x) mean(x) > 0.0001, TRUE)
root_network

#first convert to microtable
root_amp <- phyloseq2meco(root_network)
root_amp

root_env <- root_amp[["sample_table"]][, 26:46]
root_env
t1_root <- trans_env$new(dataset = root_amp, add_data = root_env)

# 'p_adjust_type = "Env"' means p adjustment is performed for each environmental variable separately.
t1_root$cal_cor(use_data = "Class", by_group = "Growing_Degree_Days", p_adjust_method = "fdr", p_adjust_type = "Env", cor_method = "pearson")

t1_root$res_cor

# default ggplot2 method with clustering
root_prop <- t1_root$plot_cor()
root_prop

root_prop[["data"]][["by_group"]] <- factor(root_prop[["data"]][["by_group"]], levels=c("600", "1400"))
root_prop



root_prop$data <- root_prop$data[
  (root_prop$data$Taxa %in% c(
    "c__Vicinamibacteria",
    "c__Verrucomicrobiae",
    "c__Thermoleophilia",
    "c__Saccharimonadia",
    "c__Polyangia",
    "c__Saccharimonadia",
    "c__Polyangia",
    "c__Planctomycetes",
    "c__Ktedonobacteria",
    "c__KD4-96",
    "c__Gemmatimonadetes",
    "c__Gammaproteobacteria",
    "c__Chloroflexia",
    "c__Blastocatellia",
    "c__Bacteroidia",
    "c__Bacilli",
    "c__Anaerolineae",
    "c__Alphaproteobacteria",
    "c__Actinobacteria",
    "c__Acidobacteriae",
    "c__Acidimicrobiia"
  )),
]
root_prop

#rename the taxaby removing the C__
root_prop$data$Taxa <- gsub("^c__", "", root_prop$data$Taxa)
root_prop


#include only the physicochemical variables that show significant correlations with class
root_prop$data <- root_prop$data[(root_prop$data$Env %in% c(
  "Ca_ppm", 
  "Cation_Exchange_Capacity",
  "Cr_ppm", 
  "Fe_ppm", 
  "K_ppm",
  "Mn_ppm",
  "Na_ppm", 
  "P_ppm"
)), 
]
root_prop

root_prop$data$Env <- factor(root_prop$data$Env, labels = c(
  "Ca_ppm" = "Calcium",
  "Cation_Exchange_Capacity" = "Cation Exchange Capacity",
  "Cr_ppm" = "Chromium",
  "Fe_ppm" = "Iron",
  "K_ppm" = "Potassium",
  "Mn_ppm" = "Manganese",
  "Na_ppm" = "Sodium",
  "P_ppm" = "Phosphorus"
))

root_prop


root_prop_final <- root_prop + theme(plot.title = element_text(family = "serif", face = "bold", size = (80)),
                                     axis.text.x=element_text(angle=45,hjust=1,vjust=1.0,size=50, colour = "black"),
                                     axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=80, colour = "black"),
                                     legend.text = element_text(size = 45, face = "bold"),
                                     strip.text = element_text(size = 70, face = "bold", angle = 0),
                                     legend.title = element_text(face="bold", size = 70),
                                     legend.key.height = unit(5, 'cm'), legend.key.width = unit(5, 'cm')) +
                                     ggtitle("Root")




root_prop_final


#save plot
ggsave("root_prop_final.png", plot =root_prop_final, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 30, height = 30, units = c("in"), device = "png")

ggsave("root_prop_final.pdf", plot = root_prop_final, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 30, height = 30, units = c("in"), device = "pdf")



#For stalk samples

stalk_bac_sub_no_bad_Filtered_final_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/stalk_bac_sub_no_bad_Filtered_final_non_normalized.rds")

stalk_norm = transform_sample_counts(stalk_bac_sub_no_bad_Filtered_final_non_normalized, function(x) x / sum(x) )
stalk_network = filter_taxa(stalk_norm, function(x) mean(x) > 0.0001, TRUE)
stalk_network

#first convert to microtable
stalk_amp <- phyloseq2meco(stalk_network)
stalk_amp

stalk_env <- stalk_amp[["sample_table"]][, 21:41]
stalk_env
t1_stalk <- trans_env$new(dataset = stalk_amp, add_data = stalk_env)

# 'p_adjust_type = "Env"' means p adjustment is performed for each environmental variable separately.
t1_stalk$cal_cor(use_data = "Class", by_group = "Growing_Degree_Days", p_adjust_method = "fdr", p_adjust_type = "Env", cor_method = "pearson")

t1_stalk$res_cor

# default ggplot2 method with clustering
stalk_prop <- t1_stalk$plot_cor()
stalk_prop

stalk_prop[["data"]][["by_group"]] <- factor(stalk_prop[["data"]][["by_group"]], levels=c("600", "1400"))
stalk_prop



stalk_prop$data <- stalk_prop$data[
  (stalk_prop$data$Taxa %in% c(
    "c__Vicinamibacteria",
    "c__Verrucomicrobiae",
    "c__TK10",
    "c__Thermoleophilia",
    "c__Subgroup_25",
    "c__Saccharimonadia",
    "c__Polyangia",
    "c__Planctomycetes",
    "c__Phycisphaerae",
    "c__KD4-96",
    "c__Gammaproteobacteria",
    "c__Deinococci",
    "c__Clostridia",
    "c__Chloroflexia",
    "c__Blastocatellia",
    "c__Bdellovibrionia",
    "c__Bacteroidia",
    "c__Bacilli",
    "c__Armatimonadia",
    "c__Alphaproteobacteria",
    "c__Actinobacteria"
  )),
]
stalk_prop

#rename the taxaby removing the C__
stalk_prop$data$Taxa <- gsub("^c__", "", stalk_prop$data$Taxa)
stalk_prop


#include only the physicochemical variables that show significant correlations with class
stalk_prop$data <- stalk_prop$data[(stalk_prop$data$Env %in% c(
  "Cr_ppm", 
  "Na_ppm",
  "Fe_ppm", 
  "Mn_ppm", 
  "pH",
  "Water_pH"
)), 
]
stalk_prop

stalk_prop$data$Env <- factor(stalk_prop$data$Env, labels = c(
  "Cr_ppm" = "Chromium",
  "Na_ppm" = "Sodium",
  "Fe_ppm" = "Iron",
  "Mn_ppm" = "Manganese",
  "pH" = "pH",
  "Water_pH" = "Water pH"
))

stalk_prop


stalk_prop_final <- stalk_prop + theme(plot.title = element_text(family = "serif", face = "bold", size = (80)),
                                       axis.text.x=element_text(angle=45,hjust=1,vjust=1.0,size=50, colour = "black"),
                                       axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=80, colour = "black"),
                                       legend.text = element_text(size = 45, face = "bold"),
                                       strip.text = element_text(size = 70, face = "bold", angle = 0),
                                       legend.title = element_text(face="bold", size = 70),
                                       legend.key.height = unit(5, 'cm'), legend.key.width = unit(5, 'cm')) +
                                       ggtitle("stalk")




stalk_prop_final


#save plot
ggsave("stalk_prop_final.png", plot =stalk_prop_final, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 30, height = 30, units = c("in"), device = "png")

ggsave("stalk_prop_final.pdf", plot = stalk_prop_final, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 30, height = 30, units = c("in"), device = "pdf")
