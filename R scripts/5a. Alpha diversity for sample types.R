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

all_samples_non_normalized_final <- readRDS(file = "C:/Aduragbemi/Microbiome/RDs/all_samples_non_normalized_final.rds")
bac.css.non_norm_root <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/root_bac_sub_no_bad_Filtered_final_non_normalized.rds")
bac.css.non_norm_stalk <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/stalk_bac_sub_no_bad_Filtered_final_non_normalized.rds")
bac.css.non_norm_soil <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/soil_bac_sub_no_bad_Filtered_final_non_normalized.rds")
bac.css.non_norm_rhizo <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/rhizo_bac_sub_no_bad_Filtered_final_non_normalized.rds")

# ===============================================================================
# Rarefy to 90% of sample depth
# Rarefaction depth chosen is the 90% of the minimum sample depth in the dataset
# ===============================================================================

all_norm_alpha_rare <- rarefy_even_depth(all_samples_non_normalized_final, sample.size=500, replace=F)
all_norm_alpha_rare

bac.css.non_norm_root_rare <- rarefy_even_depth(bac.css.non_norm_root, sample.size=1000, replace=F)
bac.css.non_norm_root_rare

bac.css.non_norm_soil_rare <- rarefy_even_depth(bac.css.non_norm_soil, sample.size=2000, replace=F)
bac.css.non_norm_soil_rare

bac.css.non_norm_stalk_rare <- rarefy_even_depth(bac.css.non_norm_stalk, sample.size=500, replace=F)
bac.css.non_norm_stalk_rare

bac.css.non_norm_rhizo_rare <- rarefy_even_depth(bac.css.non_norm_rhizo, sample.size=3000, replace=F)
bac.css.non_norm_rhizo_rare

# ===============================================================================
# Access the alpha_diversity for all sample types
# ===============================================================================

#first perform estimate richness for both observed, Shannon, faith PD's metrics
all_Sample_Type_1 <- data.frame(all_norm_alpha_rare@sam_data)
all_Sample_Type_1
all_Sample_Type_2 <- estimate_richness(all_norm_alpha_rare, measures = c("Observed", "Shannon"))
all_Sample_Type_2

# Calculate Faith's Phylogenetic Diversity (PD)
pd_all_otu <- as.data.frame(all_norm_alpha_rare@otu_table)
pd_all_tree <- all_norm_alpha_rare@phy_tree

#check if the tree is rooted
all_norm_alpha_rare@phy_tree

# it is a rooted tree
faith_pd_all <- pd(t(pd_all_otu), pd_all_tree,include.root=T)
faith_pd_all

#combine all the diversity metrics into a single data frame
all_Sample_Type_2 <- cbind(all_Sample_Type_2, all_Sample_Type_1, faith_pd_all)
all_Sample_Type_2 <- data.frame(all_Sample_Type_2)
all_Sample_Type_2

#perform shapiro normality test to check if the alpha metrics are normally distributed
shapiro.test (all_Sample_Type_2$Observed)
shapiro.test (all_Sample_Type_2$Shannon)
shapiro.test (all_Sample_Type_2$PD)

#plot the histogram to see if they are normally distributed
hist(all_Sample_Type_2$Observed, main="Observed diversity", xlab="", breaks=10)
hist(all_Sample_Type_2$Shannon, main="Shannon diversity", xlab="", breaks=10)
hist(all_Sample_Type_2$PD, main="Faith's Phylogenetic Diversity", xlab="", breaks=10)


#perform Kruskal-Wallis test  mean comparison
#for observed ASVs
kruskal.test(Observed ~ Sample_Type, data = all_Sample_Type_2)    
PT_Observed_Sample_Type_all <- dunnTest(Observed ~ Sample_Type, data = all_Sample_Type_2, method="bonferroni")     
PT2_Observed_Sample_Type_all <- PT_Observed_Sample_Type_all$res  
PT2_Observed_Sample_Type_all
PTfinal_Observed_Sample_Type_all <- cldList(comparison = PT2_Observed_Sample_Type_all$Comparison, 
                                            p.value = PT2_Observed_Sample_Type_all$P.adj, threshold  = 0.05) 
PTfinal_Observed_Sample_Type_all


#for shannon diversity
kruskal.test(Shannon ~ Sample_Type, data = all_Sample_Type_2)    
PT_Shannon_Sample_Type_all <- dunnTest(Shannon ~ Sample_Type, data = all_Sample_Type_2, method="bonferroni")      
PT2_Shannon_Sample_Type_all <- PT_Shannon_Sample_Type_all$res  
PT2_Shannon_Sample_Type_all
PTfinal_Shannon_Sample_Type_all <- cldList(comparison = PT2_Shannon_Sample_Type_all$Comparison, 
                                           p.value = PT2_Shannon_Sample_Type_all$P.adj, threshold  = 0.05) 
PTfinal_Shannon_Sample_Type_all

#for Faith's Phylogenetic Diversity (PD) 
kruskal.test(PD ~ Sample_Type, data = all_Sample_Type_2)    
PT_PD_Sample_Type_all <- dunnTest(PD ~ Sample_Type, data = all_Sample_Type_2, method="bonferroni")     
PT2_PD_Sample_Type_all <- PT_PD_Sample_Type_all$res  
PT2_PD_Sample_Type_all
PTfinal_PD_Sample_Type_all <- cldList(comparison = PT2_PD_Sample_Type_all$Comparison, 
                                      p.value = PT2_PD_Sample_Type_all$P.adj, threshold  = 0.05) 
PTfinal_PD_Sample_Type_all


#now combined all diversity metrics to data frame
all_alpha <- all_Sample_Type_2 %>%
  select(Observed, Shannon, PD, Growing_Degree_Days, Sample_Type)
all_alpha

#assign colors for the four sample types
colors <- c("#008080", "#D2691E", "#6A5ACD", "#FF1493" )

#reorder the facet grid 
all_alpha$Sample_Type <- factor(all_alpha$Sample_Type, levels=c("Soil", "Rhizosphere", "Root", "Stalk"))
all_alpha$Sample_Type

#plot for observed diversity metrics
all_final_alpha_ob <- ggplot(all_alpha, aes(x = Sample_Type, y = Observed, fill = Sample_Type))+
  geom_violin(width =0.5, aes(fill = Sample_Type), position = position_dodge(1), scale = "width") +
  scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width=0.1, fill= "white")+
  scale_y_continuous (limits = c(0, 450)) +
  ylab("Alpha Diversity Measure") +
  xlab("")+
  scale_fill_manual(values = colors) +
  theme_grey() +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", size = (15), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,size=15, colour = "black", face = "bold"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black", face = "bold"),
        legend.text = element_text(size = 15),
        legend.title = element_text(face="bold", size = 15)) +
  theme(legend.position= "") +
  ggtitle("Observed ASVs")

all_final_alpha_ob

#save as image
ggsave("all_final_alpha_ob.png", plot = all_final_alpha_ob , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("all_final_alpha_ob.pdf", plot = all_final_alpha_ob, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")


#plot for Shannon diversity metrics
all_final_alpha_sh <- ggplot(all_alpha, aes(x = Sample_Type, y = Shannon, fill = Sample_Type))+
  geom_violin(width =0.5, aes(fill = Sample_Type), position = position_dodge(1), scale = "width") +
  scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width=0.1, fill= "white")+
  scale_y_continuous (limits = c(1.5, 7)) +
  ylab("") +
  xlab("")+
  scale_fill_manual(values = colors) +
  theme_grey() +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", size = (15), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,size=15, colour = "black", face = "bold"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black", face = "bold"),
        legend.text = element_text(size = 15),
        legend.title = element_text(face="bold", size = 15)) +
  theme(legend.position="right") +
  guides(fill=guide_legend(title = "")) +
  ggtitle("Shannon Diversity")

all_final_alpha_sh

#save as image
ggsave("all_final_alpha_sh.png", plot = all_final_alpha_sh , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")


#save as pdf
ggsave("all_final_alpha_sh.pdf", plot = all_final_alpha_sh, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")


#plot for Faith's Phylogenetic Diversity (PD) metrics
all_final_alpha_pd <- ggplot(all_alpha, aes(x = Sample_Type, y = PD, fill = Sample_Type))+
  geom_violin(width =0.5, aes(fill = Sample_Type), position = position_dodge(1), scale = "width") +
  scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width=0.1, fill= "white")+
  scale_y_continuous (limits = c(0, 60)) +
  ylab("") +
  xlab("")+
  scale_fill_manual(values = colors) +
  theme_grey() +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", size = (15), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,size=15, colour = "black", face = "bold"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black", face = "bold"),
        legend.text = element_text(size = 15),
        legend.title = element_text(face="bold", size = 15)) +
  theme(legend.position="right") +
  guides(fill=guide_legend(title = "Sample Types")) +
  ggtitle("Faith's Phylogenetic Diversity")

all_final_alpha_pd

#save as image
ggsave("all_final_alpha_pd.png", plot = all_final_alpha_pd , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")


#save as pdf
ggsave("all_final_alpha_pd.pdf", plot = all_final_alpha_pd, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")



# ===============================================================================
# Access the alpha_diversity for all sample types across two growing degree days
# ===============================================================================

#Soil sample

#first perform estimate richness for both observed and Shannon
soil_gdd_1 <- data.frame(bac.css.non_norm_soil_rare@sam_data)
soil_gdd_1
soil_gdd_2 <- estimate_richness(bac.css.non_norm_soil_rare, measures = c("Observed", "Shannon"))
soil_gdd_2

#Calculate Faith's Phylogenetic Diversity (PD)
pd_soil_otu <- as.data.frame(bac.css.non_norm_soil_rare@otu_table)
pd_soil_tree <- bac.css.non_norm_soil_rare@phy_tree

#check if the tree is soiled
bac.css.non_norm_soil_rare@phy_tree

# it is a soiled tree
faith_pd_soil <- pd(t(pd_soil_otu), pd_soil_tree,include.root=T)
faith_pd_soil

#combine all the diversity metrics into a single data frame
soil_gdd_2 <- cbind(soil_gdd_2, soil_gdd_1,faith_pd_soil)
soil_gdd_2 <- data.frame(soil_gdd_2)
soil_gdd_2

#perform shapiro normality test to check if the alpha metrics are normally distributed
shapiro.test (soil_gdd_2$Observed)
shapiro.test (soil_gdd_2$Shannon)
shapiro.test (soil_gdd_2$PD)

#plot the histogram to see if they are normally distributed
hist(soil_gdd_2$Observed, main="Observed diversity", xlab="", breaks=10)
hist(soil_gdd_2$Shannon, main="Shannon diversity", xlab="", breaks=10)
hist(soil_gdd_2$PD, main="Faith's Phylogenetic Diversity", xlab="", breaks=10)

#perform Kruskal-Wallis test  mean comparison
#for observed ASvs
kruskal.test(Observed ~ Growing_Degree_Days, data = soil_gdd_2)    
PT_Observed_gdd_soil <- dunnTest(Observed ~ Growing_Degree_Days, data = soil_gdd_2, method="bonferroni")      
PT2_Observed_gdd_soil <- PT_Observed_gdd_soil$res  
PT2_Observed_gdd_soil
PTfinal_Observed_gdd_soil <- cldList(comparison = PT2_Observed_gdd_soil$Comparison, 
                                     p.value = PT2_Observed_gdd_soil$P.adj, threshold  = 0.05) 
PTfinal_Observed_gdd_soil


#for Shannon diversity
kruskal.test(Shannon ~ Growing_Degree_Days, data = soil_gdd_2)    
PT_Shannon_gdd_soil <- dunnTest(Shannon ~ Growing_Degree_Days, data = soil_gdd_2, method="bonferroni")     
PT2_Shannon_gdd_soil <- PT_Shannon_gdd_soil$res  
PT2_Shannon_gdd_soil
PTfinal_Shannon_gdd_soil <- cldList(comparison = PT2_Shannon_gdd_soil$Comparison, 
                                    p.value = PT2_Shannon_gdd_soil$P.adj, threshold  = 0.05) 
PTfinal_Shannon_gdd_soil

#for Faith's Phylogenetic Diversity 
kruskal.test(PD ~ Growing_Degree_Days, data = soil_gdd_2)    
PT_PD_gdd_soil <- dunnTest(PD ~ Growing_Degree_Days, data = soil_gdd_2, method="bonferroni")     
PT2_PD_gdd_soil <- PT_PD_gdd_soil$res  
PT2_PD_gdd_soil
PTfinal_PD_gdd_soil <- cldList(comparison = PT2_PD_gdd_soil$Comparison, 
                               p.value = PT2_PD_gdd_soil$P.adj, threshold  = 0.05) 
PTfinal_PD_gdd_soil

#combine all the diversity metrics alongside other variables into a single data frame
soil_alpha <- soil_gdd_2 %>% 
  select(Observed, Shannon, PD, Growing_Degree_Days, Sample_Type)
soil_alpha


#Rhizosphere sample

#first perform estimate richness for both observed and Shannon
rhizo_gdd_1 <- data.frame(bac.css.non_norm_rhizo_rare@sam_data)
rhizo_gdd_1
rhizo_gdd_2 <- estimate_richness(bac.css.non_norm_rhizo_rare, measures = c("Observed", "Shannon"))
rhizo_gdd_2

# Calculate Faith's Phylogenetic Diversity (PD)
pd_rhizo_otu <- as.data.frame(bac.css.non_norm_rhizo_rare@otu_table)
pd_rhizo_tree <- bac.css.non_norm_rhizo_rare@phy_tree

#check if the tree is rooted
bac.css.non_norm_rhizo_rare@phy_tree

# it is a rooted tree
faith_pd_rhizo <- pd(t(pd_rhizo_otu), pd_rhizo_tree,include.root=T)
faith_pd_rhizo

#combine all the diversity metrics into a single data frame
rhizo_gdd_2 <- cbind(rhizo_gdd_2, rhizo_gdd_1,faith_pd_rhizo)
rhizo_gdd_2 <- data.frame(rhizo_gdd_2)
rhizo_gdd_2

#perform shapiro normality test to check if the alpha metrics are normally distributed
shapiro.test (rhizo_gdd_2$Observed)
shapiro.test (rhizo_gdd_2$Shannon)
shapiro.test (rhizo_gdd_2$PD)

#plot the histogram to see if they are normally distributed
hist(rhizo_gdd_2$Observed, main="Observed diversity", xlab="", breaks=10)
hist(rhizo_gdd_2$Shannon, main="Shannon diversity", xlab="", breaks=10)
hist(rhizo_gdd_2$PD, main="Faith's Phylogenetic Diversity", xlab="", breaks=10)


#perform Kruskal-Wallis test  mean comparison
#for observed ASvs
kruskal.test(Observed ~ Growing_Degree_Days, data = rhizo_gdd_2)    
PT_Observed_gdd_rhizo <- dunnTest(Observed ~ Growing_Degree_Days, data = rhizo_gdd_2, method="bonferroni")    
PT2_Observed_gdd_rhizo <- PT_Observed_gdd_rhizo$res  
PT2_Observed_gdd_rhizo
PTfinal_Observed_gdd_rhizo <- cldList(comparison = PT2_Observed_gdd_rhizo$Comparison, 
                                      p.value = PT2_Observed_gdd_rhizo$P.adj, threshold  = 0.05) 
PTfinal_Observed_gdd_rhizo


#for Shannon diversity 
kruskal.test(Shannon ~ Growing_Degree_Days, data = rhizo_gdd_2)    
PT_Shannon_gdd_rhizo <- dunnTest(Shannon ~ Growing_Degree_Days, data = rhizo_gdd_2, method="bonferroni") 
PT2_Shannon_gdd_rhizo <- PT_Shannon_gdd_rhizo$res  
PT2_Shannon_gdd_rhizo
PTfinal_Shannon_gdd_rhizo <- cldList(comparison = PT2_Shannon_gdd_rhizo$Comparison, 
                                     p.value = PT2_Shannon_gdd_rhizo$P.adj, threshold  = 0.05) 
PTfinal_Shannon_gdd_rhizo

#for Faith's Phylogenetic Diversity 
kruskal.test(PD ~ Growing_Degree_Days, data = rhizo_gdd_2)    
PT_PD_gdd_rhizo <- dunnTest(PD ~ Growing_Degree_Days, data = rhizo_gdd_2, method="bonferroni")     
PT2_PD_gdd_rhizo <- PT_PD_gdd_rhizo$res  
PT2_PD_gdd_rhizo
PTfinal_PD_gdd_rhizo <- cldList(comparison = PT2_PD_gdd_rhizo$Comparison, 
                                p.value = PT2_PD_gdd_rhizo$P.adj, threshold  = 0.05) 
PTfinal_PD_gdd_rhizo

#combine all the diversity metrics alongside other variables into a single data frame
rhizo_alpha <- rhizo_gdd_2 %>% 
  select(Observed, Shannon, PD, Growing_Degree_Days, Sample_Type)
rhizo_alpha


#Root sample

#first perform estimate richness for both observed and Shannon
root_gdd_1 <- data.frame(bac.css.non_norm_root_rare@sam_data)
root_gdd_1
root_gdd_2 <- estimate_richness(bac.css.non_norm_root_rare, measures = c("Observed", "Shannon"))
root_gdd_2

# Calculate Faith's Phylogenetic Diversity (PD)
pd_root_otu <- as.data.frame(bac.css.non_norm_root_rare@otu_table)
pd_root_tree <- bac.css.non_norm_root_rare@phy_tree

#check if the tree is rooted
bac.css.non_norm_root_rare@phy_tree

# it is a rooted tree
faith_pd_root <- pd(t(pd_root_otu), pd_root_tree,include.root=T)
faith_pd_root

#combine all the diversity metrics into a single data frame
root_gdd_2 <- cbind(root_gdd_2, root_gdd_1,faith_pd_root)
root_gdd_2 <- data.frame(root_gdd_2)
root_gdd_2

#perform Shapiro normality test to check if the alpha metrics are normally distributed
shapiro.test (root_gdd_2$Observed)
shapiro.test (root_gdd_2$Shannon)
shapiro.test (root_gdd_2$PD)

#plot the histogram to see if they are normally distributed
hist(root_gdd_2$Observed, main="Observed diversity", xlab="", breaks=10)
hist(root_gdd_2$Shannon, main="Shannon diversity", xlab="", breaks=10)
hist(root_gdd_2$PD, main="Faith's Phylogenetic Diversity", xlab="", breaks=10)

#perform Kruskal-Wallis test  mean comparison

#for observed ASvs
kruskal.test(Observed ~ Growing_Degree_Days, data = root_gdd_2)    
PT_Observed_gdd_root <- dunnTest(Observed ~ Growing_Degree_Days, data = root_gdd_2, method="bonferroni")     
PT2_Observed_gdd_root <- PT_Observed_gdd_root$res  
PT2_Observed_gdd_root
PTfinal_Observed_gdd_root <- cldList(comparison = PT2_Observed_gdd_root$Comparison, 
                                     p.value = PT2_Observed_gdd_root$P.adj, threshold  = 0.05) 
PTfinal_Observed_gdd_root


#for shannon 
kruskal.test(Shannon ~ Growing_Degree_Days, data = root_gdd_2)    
PT_Shannon_gdd_root <- dunnTest(Shannon ~ Growing_Degree_Days, data = root_gdd_2, method="bonferroni") 
PT2_Shannon_gdd_root <- PT_Shannon_gdd_root$res  
PT2_Shannon_gdd_root
PTfinal_Shannon_gdd_root <- cldList(comparison = PT2_Shannon_gdd_root$Comparison, 
                                    p.value = PT2_Shannon_gdd_root$P.adj, threshold  = 0.05) 
PTfinal_Shannon_gdd_root

#for Faith's Phylogenetic Diversity 
kruskal.test(PD ~ Growing_Degree_Days, data = root_gdd_2)    
PT_PD_gdd_root <- dunnTest(PD ~ Growing_Degree_Days, data = root_gdd_2, method="bonferroni") 
PT2_PD_gdd_root <- PT_PD_gdd_root$res  
PT2_PD_gdd_root
PTfinal_PD_gdd_root <- cldList(comparison = PT2_PD_gdd_root$Comparison, 
                               p.value = PT2_PD_gdd_root$P.adj, threshold  = 0.05) 
PTfinal_PD_gdd_root

#combine all the diversity metrics alongside other variables into a single data frame
root_alpha <- root_gdd_2 %>% 
  select(Observed, Shannon, PD, Growing_Degree_Days, Sample_Type)
root_alpha



#for stalk sample
#first perform estimate richness for both observed and Shannon
stalk_gdd_1 <- data.frame(bac.css.non_norm_stalk_rare@sam_data)
stalk_gdd_1
stalk_gdd_2 <- estimate_richness(bac.css.non_norm_stalk_rare, measures = c("Observed", "Shannon"))
stalk_gdd_2

# Calculate Faith's Phylogenetic Diversity (PD)
pd_stalk_otu <- as.data.frame(bac.css.non_norm_stalk_rare@otu_table)
pd_stalk_tree <- bac.css.non_norm_stalk_rare@phy_tree

#check if the tree is stalked
bac.css.non_norm_stalk_rare@phy_tree

# it is a stalked tree
faith_pd_stalk <- pd(t(pd_stalk_otu), pd_stalk_tree,include.root=T)
faith_pd_stalk

#combine all the diversity metrics into a single data frame
stalk_gdd_2 <- cbind(stalk_gdd_2, stalk_gdd_1,faith_pd_stalk)
stalk_gdd_2 <- data.frame(stalk_gdd_2)
stalk_gdd_2

#perform shapiro normality test to check if the alpha metrics are normally distributed
shapiro.test (stalk_gdd_2$Observed)
shapiro.test (stalk_gdd_2$Shannon)
shapiro.test (stalk_gdd_2$PD)

#plot the histogram to see if they are normally distributed
hist(stalk_gdd_2$Observed, main="Observed diversity", xlab="", breaks=10)
hist(stalk_gdd_2$Shannon, main="Shannon diversity", xlab="", breaks=10)
hist(stalk_gdd_2$PD, main="Faith's Phylogenetic Diversity", xlab="", breaks=10)


#perform Kruskal-Wallis test  mean comparison

#for observed
kruskal.test(Observed ~ Growing_Degree_Days, data = stalk_gdd_2)    
PT_Observed_gdd_stalk <- dunnTest(Observed ~ Growing_Degree_Days, data = stalk_gdd_2, method="bonferroni")    
PT2_Observed_gdd_stalk <- PT_Observed_gdd_stalk$res  
PT2_Observed_gdd_stalk
PTfinal_Observed_gdd_stalk <- cldList(comparison = PT2_Observed_gdd_stalk$Comparison, 
                                      p.value = PT2_Observed_gdd_stalk$P.adj, threshold  = 0.05) 
PTfinal_Observed_gdd_stalk


#for shannon 
kruskal.test(Shannon ~ Growing_Degree_Days, data = stalk_gdd_2)    
PT_Shannon_gdd_stalk <- dunnTest(Shannon ~ Growing_Degree_Days, data = stalk_gdd_2, method="bonferroni")    
PT2_Shannon_gdd_stalk <- PT_Shannon_gdd_stalk$res  
PT2_Shannon_gdd_stalk
PTfinal_Shannon_gdd_stalk <- cldList(comparison = PT2_Shannon_gdd_stalk$Comparison, 
                                     p.value = PT2_Shannon_gdd_stalk$P.adj, threshold  = 0.05) 
PTfinal_Shannon_gdd_stalk

#for Faith's Phylogenetic Diversity 
kruskal.test(PD ~ Growing_Degree_Days, data = stalk_gdd_2)    
PT_PD_gdd_stalk <- dunnTest(PD ~ Growing_Degree_Days, data = stalk_gdd_2, method="bonferroni")    
PT2_PD_gdd_stalk <- PT_PD_gdd_stalk$res  
PT2_PD_gdd_stalk
PTfinal_PD_gdd_stalk <- cldList(comparison = PT2_PD_gdd_stalk$Comparison, 
                                p.value = PT2_PD_gdd_stalk$P.adj, threshold  = 0.05) 
PTfinal_PD_gdd_stalk

#combine all the diversity metrics alongside other variables into a single data frame
stalk_alpha <- stalk_gdd_2 %>% 
  select(Observed, Shannon, PD, Growing_Degree_Days, Sample_Type)
stalk_alpha


#merge all together
combined_alpha <- bind_rows(soil_alpha, rhizo_alpha, root_alpha, stalk_alpha)
combined_alpha

phylum_colors <- c("#7FCDBB", "#006D6F")


#reorder the facet grid 
combined_alpha$Sample_Type <- factor(combined_alpha$Sample_Type, levels=c("Soil", "Rhizosphere", "Root", "Stalk"))
combined_alpha$Growing_Degree_Days <- factor(combined_alpha$Growing_Degree_Days, levels=c("600", "1400"))


#plot for observed diversity metrics
combined_final_alpha_ob <- ggplot(combined_alpha, aes(x = Growing_Degree_Days, y = Observed, fill = Growing_Degree_Days))+
  facet_grid(cols = vars(Sample_Type))+
  geom_violin(width =0.5, aes(fill = Growing_Degree_Days), position = position_dodge(1), scale = "width") +
  scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width=0.1, fill= "white", outlier.colour ="NA")+
  scale_y_continuous (limits = c(0, 900)) +
  ylab("Alpha Diversity Measure") +
  xlab("")+
  scale_fill_manual(values = phylum_colors) +
  theme_grey() +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", size = (15), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,size=15, colour = "black"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black"),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 15)) +
  theme(legend.position="right") +
  guides(fill = guide_legend("")) +
  ggtitle("Observed ASVs")

combined_final_alpha_ob

#save as image
ggsave("combined_final_alpha_ob.png", plot = combined_final_alpha_ob , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("combined_final_alpha_ob.pdf", plot = combined_final_alpha_ob, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")


#plot for Shannon diversity metrics
#reorder the facet grid 
combined_alpha$Sample_Type <- factor(combined_alpha$Sample_Type, levels=c("Soil", "Rhizosphere", "Root", "Stalk"))
combined_alpha$Growing_Degree_Days <- factor(combined_alpha$Growing_Degree_Days, levels=c("600", "1400"))

combined_final_alpha_sh <- ggplot(combined_alpha, aes(x = Growing_Degree_Days, y = Shannon, fill = Growing_Degree_Days))+
  facet_grid(cols = vars(Sample_Type))+
  geom_violin(width =0.5, aes(fill = Growing_Degree_Days), position = position_dodge(1), scale = "width") +
  scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width=0.1, fill= "white")+
  ylab("") +
  xlab("Growing Degree Days")+
  scale_y_continuous (limits = c(1.5, 7))  +
  scale_fill_manual(values = phylum_colors) +
  theme_grey() +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", size = (15), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,size=15, colour = "black"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black"),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 15)) +
  theme(legend.position="right") +
  guides(fill = guide_legend("")) +
  ggtitle("Shannon Diversity")

combined_final_alpha_sh

#save as image
ggsave("combined_final_alpha_sh.png", plot = combined_final_alpha_sh , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("combined_final_alpha_sh.pdf", plot = combined_final_alpha_sh, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")


#plot for Faith's Phylogenetic Diversity metrics
#reorder the facet grid 
combined_alpha$Sample_Type <- factor(combined_alpha$Sample_Type, levels=c("Soil", "Rhizosphere", "Root", "Stalk"))
combined_alpha$Growing_Degree_Days <- factor(combined_alpha$Growing_Degree_Days, levels=c("600", "1400"))

combined_final_alpha_pd <- ggplot(combined_alpha, aes(x = Growing_Degree_Days, y = PD, fill = Growing_Degree_Days))+
  facet_grid(cols = vars(Sample_Type))+
  geom_violin(width =0.5, aes(fill = Growing_Degree_Days), position = position_dodge(1), scale = "width") +
  scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width=0.1, fill= "white")+
  ylab("") +
  xlab("Growing Degree Days")+
  scale_y_continuous (limits = c(0, 75))  +
  scale_fill_manual(values = phylum_colors) +
  theme_grey() +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", size = (15), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,size=15, colour = "black"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black"),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 15)) +
  theme(legend.position="right") +
  guides(fill = guide_legend("Growing Degree Days")) +
  ggtitle("Faith's Phylogenetic Diversity")

combined_final_alpha_pd


#save as image
ggsave("combined_final_alpha_pd.png", plot = combined_final_alpha_pd , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("combined_final_alpha_pd.pdf", plot = combined_final_alpha_pd, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")




