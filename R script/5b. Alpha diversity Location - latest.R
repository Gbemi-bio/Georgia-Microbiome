#Load Libraries
library(FSA) 
library(rcompanion)
library(agricolae)
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(picante)


# Set seed for reproducible result
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

# ==========================================================================
# alpha-diversity measures of bacterial communities across various locations
# ==========================================================================

#Soil sample
#first perform estimate richness for both observed and Shannon
soil_location_1 <- data.frame(bac.css.non_norm_soil_rare@sam_data)
soil_location_1
soil_location_2 <- estimate_richness(bac.css.non_norm_soil_rare, measures = c("Observed", "Shannon"))
soil_location_2


#Calculate Faith's Phylogenetic Diversity (PD)
pd_soil_otu_location <- as.data.frame(bac.css.non_norm_soil_rare@otu_table)
pd_soil_tree_location <- bac.css.non_norm_soil_rare@phy_tree

#check if the tree is rooted
bac.css.non_norm_soil_rare@phy_tree

#it is a rooted tree
faith_pd_soil_location <- pd(t(pd_soil_otu_location), pd_soil_tree_location,include.root=T)
faith_pd_soil_location

#combine all the diversity metrics into a single data frame
soil_location_2 <- cbind(soil_location_2, soil_location_1, faith_pd_soil_location)
soil_location_2 <- data.frame(soil_location_2)
soil_location_2


#perform shapiro normality test to check if the alpha metrics are normally distributed
shapiro.test (soil_location_2$Observed)
shapiro.test (soil_location_2$Shannon)
shapiro.test (soil_location_2$PD)



#plot the histogram to see if they are normally distributed
hist(soil_location_2$Observed, main="Observed diversity", xlab="", breaks=10)
hist(soil_location_2$Shannon, main="Shannon diversity", xlab="", breaks=10)
hist(soil_location_2$PD, main="Faith's Phylogenetic Diversity", xlab="", breaks=10)

#perform Kruskal-Wallis test  mean comparison
#for observed ASVs
kruskal.test(Observed ~ FieldID, data = soil_location_2)    
PT_Observed_location <- dunnTest(Observed ~ FieldID, data = soil_location_2, method="bonferroni") 
PT2_Observed_location <- PT_Observed_location$res  
PT2_Observed_location

#for Shannon 
kruskal.test(Shannon ~ FieldID, data = soil_location_2)    
PT_Shannon_location <- dunnTest(Shannon ~ FieldID, data = soil_location_2, method="bonferroni")  
PT2_Shannon_location <- PT_Shannon_location$res  
PT2_Shannon_location

#for Faith's Phylogenetic Diversity (PD) 
kruskal.test(PD ~ FieldID, data = soil_location_2)    
PT_PD_location <- dunnTest(PD ~ FieldID, data = soil_location_2, method="bonferroni")     
PT2_PD_location <- PT_PD_location$res  
PT2_PD_location



#now combined all diversity metrics to data frame
soil_location_alpha <- soil_location_2 %>%
  select(Observed, Shannon, PD, Growing_Degree_Days, Sample_Type, FieldID)
soil_location_alpha



#Rhizosphere sample
#first perform estimate richness for both observed and Shannon
rhizo_location_1 <- data.frame(bac.css.non_norm_rhizo_rare@sam_data)
rhizo_location_1
rhizo_location_2 <- estimate_richness(bac.css.non_norm_rhizo_rare, measures = c("Observed", "Shannon"))
rhizo_location_2


#Calculate Faith's Phylogenetic Diversity (PD)
pd_rhizo_otu_location <- as.data.frame(bac.css.non_norm_rhizo_rare@otu_table)
pd_rhizo_tree_location <- bac.css.non_norm_rhizo_rare@phy_tree

#check if the tree is rooted
bac.css.non_norm_rhizo_rare@phy_tree

#it is a rooted tree
faith_pd_rhizo_location <- pd(t(pd_rhizo_otu_location), pd_rhizo_tree_location,include.root=T)
faith_pd_rhizo_location

#combine all the diversity metrics into a single data frame
rhizo_location_2 <- cbind(rhizo_location_2, rhizo_location_1, faith_pd_rhizo_location)
rhizo_location_2 <- data.frame(rhizo_location_2)
rhizo_location_2


#perform shapiro normality test to check if the alpha metrics are normally distributed
shapiro.test (rhizo_location_2$Observed)
shapiro.test (rhizo_location_2$Shannon)
shapiro.test (rhizo_location_2$PD)



#plot the histogram to see if they are normally distributed
hist(rhizo_location_2$Observed, main="Observed diversity", xlab="", breaks=10)
hist(rhizo_location_2$Shannon, main="Shannon diversity", xlab="", breaks=10)
hist(rhizo_location_2$PD, main="Faith's Phylogenetic Diversity", xlab="", breaks=10)

#perform Kruskal-Wallis test  mean comparison
#for observed ASVs
kruskal.test(Observed ~ FieldID, data = rhizo_location_2)    
PT_Observed_location <- dunnTest(Observed ~ FieldID, data = rhizo_location_2, method="bonferroni") 
PT2_Observed_location <- PT_Observed_location$res  
PT2_Observed_location

#for Shannon 
kruskal.test(Shannon ~ FieldID, data = rhizo_location_2)    
PT_Shannon_location <- dunnTest(Shannon ~ FieldID, data = rhizo_location_2, method="bonferroni")  
PT2_Shannon_location <- PT_Shannon_location$res  
PT2_Shannon_location

#for Faith's Phylogenetic Diversity (PD) 
kruskal.test(PD ~ FieldID, data = rhizo_location_2)    
PT_PD_location <- dunnTest(PD ~ FieldID, data = rhizo_location_2, method="bonferroni")     
PT2_PD_location <- PT_PD_location$res  
PT2_PD_location



#now combined all diversity metrics to data frame
rhizo_location_alpha <- rhizo_location_2 %>%
  select(Observed, Shannon, PD, Growing_Degree_Days, Sample_Type, FieldID)
rhizo_location_alpha


#Root sample
#first perform estimate richness for both observed and Shannon
root_location_1 <- data.frame(bac.css.non_norm_root_rare@sam_data)
root_location_1
root_location_2 <- estimate_richness(bac.css.non_norm_root_rare, measures = c("Observed", "Shannon"))
root_location_2


#Calculate Faith's Phylogenetic Diversity (PD)
pd_root_otu_location <- as.data.frame(bac.css.non_norm_root_rare@otu_table)
pd_root_tree_location <- bac.css.non_norm_root_rare@phy_tree

#check if the tree is rooted
bac.css.non_norm_root_rare@phy_tree

#it is a rooted tree
faith_pd_root_location <- pd(t(pd_root_otu_location), pd_root_tree_location,include.root=T)
faith_pd_root_location

#combine all the diversity metrics into a single data frame
root_location_2 <- cbind(root_location_2, root_location_1, faith_pd_root_location)
root_location_2 <- data.frame(root_location_2)
root_location_2


#perform shapiro normality test to check if the alpha metrics are normally distributed
shapiro.test (root_location_2$Observed)
shapiro.test (root_location_2$Shannon)
shapiro.test (root_location_2$PD)



#plot the histogram to see if they are normally distributed
hist(root_location_2$Observed, main="Observed diversity", xlab="", breaks=10)
hist(root_location_2$Shannon, main="Shannon diversity", xlab="", breaks=10)
hist(root_location_2$PD, main="Faith's Phylogenetic Diversity", xlab="", breaks=10)

#perform Kruskal-Wallis test  mean comparison
#for observed ASVs
kruskal.test(Observed ~ FieldID, data = root_location_2)    
PT_Observed_location <- dunnTest(Observed ~ FieldID, data = root_location_2, method="bonferroni") 
PT2_Observed_location <- PT_Observed_location$res  
PT2_Observed_location

#for Shannon 
kruskal.test(Shannon ~ FieldID, data = root_location_2)    
PT_Shannon_location <- dunnTest(Shannon ~ FieldID, data = root_location_2, method="bonferroni")  
PT2_Shannon_location <- PT_Shannon_location$res  
PT2_Shannon_location

#for Faith's Phylogenetic Diversity (PD) 
kruskal.test(PD ~ FieldID, data = root_location_2)    
PT_PD_location <- dunnTest(PD ~ FieldID, data = root_location_2, method="bonferroni")     
PT2_PD_location <- PT_PD_location$res  
PT2_PD_location


#now combined all diversity metrics to data frame
root_location_alpha <- root_location_2 %>%
  select(Observed, Shannon, PD, Growing_Degree_Days, Sample_Type, FieldID)
root_location_alpha




#Stalk sample
#first perform estimate richness for both observed and Shannon
stalk_location_1 <- data.frame(bac.css.non_norm_stalk_rare@sam_data)
stalk_location_1
stalk_location_2 <- estimate_richness(bac.css.non_norm_stalk_rare, measures = c("Observed", "Shannon"))
stalk_location_2


#Calculate Faith's Phylogenetic Diversity (PD)
pd_stalk_otu_location <- as.data.frame(bac.css.non_norm_stalk_rare@otu_table)
pd_stalk_tree_location <- bac.css.non_norm_stalk_rare@phy_tree

#check if the tree is rooted
bac.css.non_norm_stalk_rare@phy_tree

#it is a rooted tree
faith_pd_stalk_location <- pd(t(pd_stalk_otu_location), pd_stalk_tree_location,include.root=T)
faith_pd_stalk_location

#combine all the diversity metrics into a single data frame
stalk_location_2 <- cbind(stalk_location_2, stalk_location_1, faith_pd_stalk_location)
stalk_location_2 <- data.frame(stalk_location_2)
stalk_location_2


#perform shapiro normality test to check if the alpha metrics are normally distributed
shapiro.test (stalk_location_2$Observed)
shapiro.test (stalk_location_2$Shannon)
shapiro.test (stalk_location_2$PD)



#plot the histogram to see if they are normally distributed
hist(stalk_location_2$Observed, main="Observed diversity", xlab="", breaks=10)
hist(stalk_location_2$Shannon, main="Shannon diversity", xlab="", breaks=10)
hist(stalk_location_2$PD, main="Faith's Phylogenetic Diversity", xlab="", breaks=10)

#perform Kruskal-Wallis test  mean comparison
#for observed ASVs
kruskal.test(Observed ~ FieldID, data = stalk_location_2)    
PT_Observed_location <- dunnTest(Observed ~ FieldID, data = stalk_location_2, method="bonferroni") 
PT2_Observed_location <- PT_Observed_location$res  
PT2_Observed_location

#for Shannon 
kruskal.test(Shannon ~ FieldID, data = stalk_location_2)    
PT_Shannon_location <- dunnTest(Shannon ~ FieldID, data = stalk_location_2, method="bonferroni")  
PT2_Shannon_location <- PT_Shannon_location$res  
PT2_Shannon_location

#for Faith's Phylogenetic Diversity (PD) 
kruskal.test(PD ~ FieldID, data = stalk_location_2)    
PT_PD_location <- dunnTest(PD ~ FieldID, data = stalk_location_2, method="bonferroni")     
PT2_PD_location <- PT_PD_location$res  
PT2_PD_location



#now combined all diversity metrics to data frame
stalk_location_alpha <- stalk_location_2 %>%
  select(Observed, Shannon, PD, Growing_Degree_Days, Sample_Type, FieldID)
stalk_location_alpha



#=============================================================
# Function to create publication-quality alpha diversity plots
#=============================================================

plot_alpha <- function(data,
                       response,
                       ylab,
                       title,
                       ylim,
                       save_name,
                       image_path,
                       pdf_path){
  

  # Summary data for vertical lines
  summary_df <- data %>%
    group_by(FieldID) %>%
    summarise(
      ymin = min(.data[[response]], na.rm = TRUE),
      ymax = max(.data[[response]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # Kruskal-Wallis p-value
  kw <- kruskal.test(
    as.formula(paste(response, "~ FieldID")),
    data = data
  )
  
  p_label <- paste0(
    "p-value = ",
    signif(kw$p.value,3)
  )
  
  # Plot
  p <- ggplot() +
    
    geom_segment(
      data = summary_df,
      aes(
        x = FieldID,
        xend = FieldID,
        y = ymin,
        yend = ymax
      ),
      colour = "grey30",
      linewidth = 0.8
    ) +
    
    geom_jitter(
      data = data,
      aes(
        x = FieldID,
        y = .data[[response]]
      ),
      width = 0.18,
      height = 0,
      size = 8,
      shape = 16,
      colour = "black",
      alpha = 0.85
    ) +
    
    scale_y_continuous(
      limits = ylim,
      expand = expansion(mult = c(0.02,0.05))
    ) +
    
    labs(
      title = title,
      x = "",
      y = ylab
    ) +
    
    annotate(
      "text",
      x = 2,
      y = ylim[2]*0.97,
      label = p_label,
      hjust = 0,
      size = 10,
      fontface = "bold"
    ) +
    
    theme_bw() +
    
    theme(
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      plot.title = element_text(
        family="serif",
        face="bold",
        size=18,
        hjust=0.5
      ),
      
      axis.title = element_text(
        family="serif",
        face="bold",
        size=16
      ),
      
      axis.text.x = element_text(
        angle=60,
        hjust=1,
        size=20,
        colour="black"
      ),
      
      axis.text.y = element_text(
        size=20,
        colour="black"
      ),
      
      axis.line = element_line(colour="black"),
      
      legend.position="none"
    )
  
  print(p)
  
  ggsave(
    paste0(save_name,".png"),
    plot=p,
    path=image_path,
    dpi=700,
    width=15,
    height=10,
    units="in"
  )
  
  ggsave(
    paste0(save_name,".pdf"),
    plot=p,
    path=pdf_path,
    dpi=700,
    width=15,
    height=10,
    units="in"
  )
  
  return(p)
  
}


#For Soil sample
soil_location_alpha_ob <-
  plot_alpha(
    data = soil_location_alpha,
    response = "Observed",
    ylab = "Observed ASVs",
    title = "Soil - Observed ASVs",
    ylim = c(0,1000),
    save_name = "soil_location_alpha_ob",
    image_path = "C:/Aduragbemi/Manuscript/Review 7/IMAGES7",
    pdf_path = "C:/Aduragbemi/Manuscript/Review 7/PDF7"
  )

soil_location_alpha_ob

soil_location_alpha_sh <-
  plot_alpha(
    data = soil_location_alpha,
    response = "Shannon",
    ylab = "Shannon Diversity",
    title = "Shannon Diversity",
    ylim = c(2,8),
    save_name = "soil_location_alpha_sh",
    image_path = "C:/Aduragbemi/Manuscript/Review 7/IMAGES7",
    pdf_path = "C:/Aduragbemi/Manuscript/Review 7/PDF7"
  )
soil_location_alpha_sh

soil_location_alpha_pd <-
  plot_alpha(
    data = soil_location_alpha,
    response = "PD",
    ylab = "Faith's PD",
    title = "Faith's Phylogenetic Diversity",
    ylim = c(0,80),
    save_name = "soil_location_alpha_pd",
    image_path = "C:/Aduragbemi/Manuscript/Review 7/IMAGES7",
    pdf_path = "C:/Aduragbemi/Manuscript/Review 7/PDF7"
  )


soil_location_alpha_pd




#For Rhizosphere sample
rhizo_location_alpha_ob <-
  plot_alpha(
    data = rhizo_location_alpha,
    response = "Observed",
    ylab = "Observed ASVs",
    title = "rhizo - Observed ASVs",
    ylim = c(0,1000),
    save_name = "rhizo_location_alpha_ob",
    image_path = "C:/Aduragbemi/Manuscript/Review 7/IMAGES7",
    pdf_path = "C:/Aduragbemi/Manuscript/Review 7/PDF7"
  )

rhizo_location_alpha_ob

rhizo_location_alpha_sh <-
  plot_alpha(
    data = rhizo_location_alpha,
    response = "Shannon",
    ylab = "Shannon Diversity",
    title = "Shannon Diversity",
    ylim = c(3,7),
    save_name = "rhizo_location_alpha_sh",
    image_path = "C:/Aduragbemi/Manuscript/Review 7/IMAGES7",
    pdf_path = "C:/Aduragbemi/Manuscript/Review 7/PDF7"
  )
rhizo_location_alpha_sh

rhizo_location_alpha_pd <-
  plot_alpha(
    data = rhizo_location_alpha,
    response = "PD",
    ylab = "Faith's PD",
    title = "Faith's Phylogenetic Diversity",
    ylim = c(0,80),
    save_name = "rhizo_location_alpha_pd",
    image_path = "C:/Aduragbemi/Manuscript/Review 7/IMAGES7",
    pdf_path = "C:/Aduragbemi/Manuscript/Review 7/PDF7"
  )


rhizo_location_alpha_pd



#For Root sample
root_location_alpha_ob <-
  plot_alpha(
    data = root_location_alpha,
    response = "Observed",
    ylab = "Observed ASVs",
    title = "root - Observed ASVs",
    ylim = c(0,500),
    save_name = "root_location_alpha_ob",
    image_path = "C:/Aduragbemi/Manuscript/Review 7/IMAGES7",
    pdf_path = "C:/Aduragbemi/Manuscript/Review 7/PDF7"
  )

root_location_alpha_ob

root_location_alpha_sh <-
  plot_alpha(
    data = root_location_alpha,
    response = "Shannon",
    ylab = "Shannon Diversity",
    title = "Shannon Diversity",
    ylim = c(3,7),
    save_name = "root_location_alpha_sh",
    image_path = "C:/Aduragbemi/Manuscript/Review 7/IMAGES7",
    pdf_path = "C:/Aduragbemi/Manuscript/Review 7/PDF7"
  )
root_location_alpha_sh

root_location_alpha_pd <-
  plot_alpha(
    data = root_location_alpha,
    response = "PD",
    ylab = "Faith's PD",
    title = "Faith's Phylogenetic Diversity",
    ylim = c(0,50),
    save_name = "root_location_alpha_pd",
    image_path = "C:/Aduragbemi/Manuscript/Review 7/IMAGES7",
    pdf_path = "C:/Aduragbemi/Manuscript/Review 7/PDF7"
  )


root_location_alpha_pd



#For Stalk sample
stalk_location_alpha_ob <-
  plot_alpha(
    data = stalk_location_alpha,
    response = "Observed",
    ylab = "Observed ASVs",
    title = "Stalk - Observed ASVs",
    ylim = c(0,150),
    save_name = "stalk_location_alpha_ob",
    image_path = "C:/Aduragbemi/Manuscript/Review 7/IMAGES7",
    pdf_path = "C:/Aduragbemi/Manuscript/Review 7/PDF7"
  )

stalk_location_alpha_ob

stalk_location_alpha_sh <-
  plot_alpha(
    data = stalk_location_alpha,
    response = "Shannon",
    ylab = "Shannon Diversity",
    title = "Shannon Diversity",
    ylim = c(1,5),
    save_name = "stalk_location_alpha_sh",
    image_path = "C:/Aduragbemi/Manuscript/Review 7/IMAGES7",
    pdf_path = "C:/Aduragbemi/Manuscript/Review 7/PDF7"
  )
stalk_location_alpha_sh


stalk_location_alpha_pd <-
  plot_alpha(
    data = stalk_location_alpha,
    response = "PD",
    ylab = "Faith's PD",
    title = "Faith's Phylogenetic Diversity",
    ylim = c(0,20),
    save_name = "stalk_location_alpha_pd",
    image_path = "C:/Aduragbemi/Manuscript/Review 7/IMAGES7",
    pdf_path = "C:/Aduragbemi/Manuscript/Review 7/PDF7"
  )

stalk_location_alpha_pd

