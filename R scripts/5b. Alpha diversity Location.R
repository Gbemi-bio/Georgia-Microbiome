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


#Make final plot for all diversity metrics
#plot for observed diversity metrics
soil_location_alpha_ob <- ggplot(soil_location_alpha, aes(x = FieldID, y = Observed)) +
  geom_violin(width =0.7, position = position_dodge(1), scale = "width") + 
  scale_fill_brewer(palette="Dark2") + 
  geom_boxplot(width= 0.1, fill= "white")+
  scale_y_continuous (limits = c(0, 1000)) +
  ylab("Alpha Diversity Measure") +
  xlab("")+
  scale_fill_manual(values = colors) +
  theme_grey() +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", size = (15), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=15, colour = "black"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black"),
        legend.text = element_text(size = 15),
        legend.title = element_text(face="bold", size = 15)) +
  theme(legend.position="") +
  ggtitle("soil - Observed ASVs")

soil_location_alpha_ob

#save as image
ggsave("soil_location_alpha_ob.png", plot = soil_location_alpha_ob, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("soil_location_alpha_ob.pdf", plot = soil_location_alpha_ob, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")

#plot for Shannon diversity metrics
soil_location_alpha_sh <- ggplot(soil_location_alpha, aes(x = FieldID, y = Shannon))+
  geom_violin(width =0.7, position = position_dodge(1), scale = "width") +
  scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width=0.1, fill= "white")+
  scale_y_continuous (limits = c(2, 8)) +
  ylab("") +
  xlab("") +
  theme_grey() +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", size = (15), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=15, colour = "black"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black"),
        legend.text = element_text(size = 15),
        legend.title = element_text(face="bold", size = 15)) +
  theme(legend.position="") +
  ggtitle("Shannon Diversity")

soil_location_alpha_sh

#save as image
ggsave("soil_location_alpha_sh.png", plot = soil_location_alpha_sh , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("soil_location_alpha_sh.pdf", plot = soil_location_alpha_sh, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")



#plot for Faith's Phylogenetic Diversity (PD) metrics
soil_location_alpha_pd <- ggplot(soil_location_alpha, aes(x = FieldID, y = PD))+
  geom_violin(width =0.7, position = position_dodge(1), scale = "width") +
  scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width=0.1, fill= "white")+
  scale_y_continuous (limits = c(0, 80)) +
  ylab("") +
  xlab("") +
  theme_grey() +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", size = (15), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=15, colour = "black"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black"),
        legend.text = element_text(size = 15),
        legend.title = element_text(face="bold", size = 15)) +
  theme(legend.position="") +
  ggtitle("Faith's Phylogenetic Diversity")

soil_location_alpha_pd

#save as image
ggsave("soil_location_alpha_pd.png", plot = soil_location_alpha_pd, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("soil_location_alpha_pd.pdf", plot = soil_location_alpha_pd, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")




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


#Make final plot for all diversity metrics
#plot for observed diversity metrics
rhizo_location_alpha_ob <- ggplot(rhizo_location_alpha, aes(x = FieldID, y = Observed)) +
  geom_violin(width =0.7, position = position_dodge(1), scale = "width") + 
  scale_fill_brewer(palette="Dark2") + 
  geom_boxplot(width= 0.1, fill= "white")+
  scale_y_continuous (limits = c(0, 1000)) +
  ylab("Alpha Diversity Measure") +
  xlab("")+
  scale_fill_manual(values = colors) +
  theme_grey() +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", size = (15), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=15, colour = "black"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black"),
        legend.text = element_text(size = 15),
        legend.title = element_text(face="bold", size = 15)) +
  theme(legend.position="") +
  ggtitle("rhizo - Observed ASVs")

rhizo_location_alpha_ob

#save as image
ggsave("rhizo_location_alpha_ob.png", plot = rhizo_location_alpha_ob, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("rhizo_location_alpha_ob.pdf", plot = rhizo_location_alpha_ob, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")

#plot for Shannon diversity metrics
rhizo_location_alpha_sh <- ggplot(rhizo_location_alpha, aes(x = FieldID, y = Shannon))+
  geom_violin(width =0.7, position = position_dodge(1), scale = "width") +
  scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width=0.1, fill= "white")+
  scale_y_continuous (limits = c(2, 8)) +
  ylab("") +
  xlab("") +
  theme_grey() +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", size = (15), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=15, colour = "black"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black"),
        legend.text = element_text(size = 15),
        legend.title = element_text(face="bold", size = 15)) +
  theme(legend.position="") +
  ggtitle("Shannon Diversity")

rhizo_location_alpha_sh

#save as image
ggsave("rhizo_location_alpha_sh.png", plot = rhizo_location_alpha_sh , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("rhizo_location_alpha_sh.pdf", plot = rhizo_location_alpha_sh, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")



#plot for Faith's Phylogenetic Diversity (PD) metrics
rhizo_location_alpha_pd <- ggplot(rhizo_location_alpha, aes(x = FieldID, y = PD))+
  geom_violin(width =0.7, position = position_dodge(1), scale = "width") +
  scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width=0.1, fill= "white")+
  scale_y_continuous (limits = c(0, 80)) +
  ylab("") +
  xlab("") +
  theme_grey() +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", size = (15), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=15, colour = "black"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black"),
        legend.text = element_text(size = 15),
        legend.title = element_text(face="bold", size = 15)) +
  theme(legend.position="") +
  ggtitle("Faith's Phylogenetic Diversity")

rhizo_location_alpha_pd

#save as image
ggsave("rhizo_location_alpha_pd.png", plot = rhizo_location_alpha_pd, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("rhizo_location_alpha_pd.pdf", plot = rhizo_location_alpha_pd, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")




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


#Make final plot for all diversity metrics
#plot for observed diversity metrics
root_location_alpha_ob <- ggplot(root_location_alpha, aes(x = FieldID, y = Observed)) +
  geom_violin(width =0.7, position = position_dodge(1), scale = "width") + 
  scale_fill_brewer(palette="Dark2") + 
  geom_boxplot(width= 0.1, fill= "white")+
  scale_y_continuous (limits = c(0, 600)) +
  ylab("Alpha Diversity Measure") +
  xlab("")+
  scale_fill_manual(values = colors) +
  theme_grey() +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", size = (15), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=15, colour = "black"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black"),
        legend.text = element_text(size = 15),
        legend.title = element_text(face="bold", size = 15)) +
  theme(legend.position="") +
  ggtitle("Root - Observed ASVs")

root_location_alpha_ob

#save as image
ggsave("root_location_alpha_ob.png", plot = root_location_alpha_ob, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("root_location_alpha_ob.pdf", plot = root_location_alpha_ob, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")

#plot for Shannon diversity metrics
root_location_alpha_sh <- ggplot(root_location_alpha, aes(x = FieldID, y = Shannon))+
  geom_violin(width =0.7, position = position_dodge(1), scale = "width") +
  scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width=0.1, fill= "white")+
  scale_y_continuous (limits = c(2, 7)) +
  ylab("") +
  xlab("") +
  theme_grey() +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", size = (15), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=15, colour = "black"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black"),
        legend.text = element_text(size = 15),
        legend.title = element_text(face="bold", size = 15)) +
  theme(legend.position="") +
  ggtitle("Shannon Diversity")

root_location_alpha_sh

#save as image
ggsave("root_location_alpha_sh.png", plot = root_location_alpha_sh , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("root_location_alpha_sh.pdf", plot = root_location_alpha_sh, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")



#plot for Faith's Phylogenetic Diversity (PD) metrics
root_location_alpha_pd <- ggplot(root_location_alpha, aes(x = FieldID, y = PD))+
  geom_violin(width =0.7, position = position_dodge(1), scale = "width") +
  scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width=0.1, fill= "white")+
  scale_y_continuous (limits = c(0, 50)) +
  ylab("") +
  xlab("") +
  theme_grey() +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", size = (15), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=15, colour = "black"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black"),
        legend.text = element_text(size = 15),
        legend.title = element_text(face="bold", size = 15)) +
  theme(legend.position="") +
  ggtitle("Faith's Phylogenetic Diversity")

root_location_alpha_pd

#save as image
ggsave("root_location_alpha_pd.png", plot = root_location_alpha_pd, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("root_location_alpha_pd.pdf", plot = root_location_alpha_pd, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")


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


#Make final plot for all diversity metrics
#plot for observed diversity metrics
stalk_location_alpha_ob <- ggplot(stalk_location_alpha, aes(x = FieldID, y = Observed)) +
  geom_violin(width =0.7, position = position_dodge(1), scale = "width") + 
  scale_fill_brewer(palette="Dark2") + 
  geom_boxplot(width= 0.1, fill= "white")+
  scale_y_continuous (limits = c(0, 150)) +
  ylab("Alpha Diversity Measure") +
  xlab("")+
  scale_fill_manual(values = colors) +
  theme_grey() +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", size = (15), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=15, colour = "black"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black"),
        legend.text = element_text(size = 15),
        legend.title = element_text(face="bold", size = 15)) +
  theme(legend.position="") +
  ggtitle("stalk - Observed ASVs")

stalk_location_alpha_ob

#save as image
ggsave("stalk_location_alpha_ob.png", plot = stalk_location_alpha_ob, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("stalk_location_alpha_ob.pdf", plot = stalk_location_alpha_ob, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")

#plot for Shannon diversity metrics
stalk_location_alpha_sh <- ggplot(stalk_location_alpha, aes(x = FieldID, y = Shannon))+
  geom_violin(width =0.7, position = position_dodge(1), scale = "width") +
  scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width=0.1, fill= "white")+
  scale_y_continuous (limits = c(2, 6)) +
  ylab("") +
  xlab("") +
  theme_grey() +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", size = (15), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=15, colour = "black"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black"),
        legend.text = element_text(size = 15),
        legend.title = element_text(face="bold", size = 15)) +
  theme(legend.position="") +
  ggtitle("Shannon Diversity")

stalk_location_alpha_sh

#save as image
ggsave("stalk_location_alpha_sh.png", plot = stalk_location_alpha_sh , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("stalk_location_alpha_sh.pdf", plot = stalk_location_alpha_sh, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")



#plot for Faith's Phylogenetic Diversity (PD) metrics
stalk_location_alpha_pd <- ggplot(stalk_location_alpha, aes(x = FieldID, y = PD))+
  geom_violin(width =0.7, position = position_dodge(1), scale = "width") +
  scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width=0.1, fill= "white")+
  scale_y_continuous (limits = c(0, 20)) +
  ylab("") +
  xlab("") +
  theme_grey() +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", size = (15), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=15, colour = "black"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black"),
        legend.text = element_text(size = 15),
        legend.title = element_text(face="bold", size = 15)) +
  theme(legend.position="") +
  ggtitle("Faith's Phylogenetic Diversity")

stalk_location_alpha_pd

#save as image
ggsave("stalk_location_alpha_pd.png", plot = stalk_location_alpha_pd, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("stalk_location_alpha_pd.pdf", plot = stalk_location_alpha_pd, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")

