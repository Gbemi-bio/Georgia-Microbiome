#Load Libraries
library(dada2)
library(phyloseq)
library(Biostrings)
library(DECIPHER)
library(phangorn)
library(dplyr)
library(ggplot2)


# ================================
# read in the non-normalized data
# ================================
root_bac_sub_no_bad_Filtered_final_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/root_bac_sub_no_bad_Filtered_final_non_normalized.rds")
stalk_bac_sub_no_bad_Filtered_final_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/stalk_bac_sub_no_bad_Filtered_final_non_normalized.rds")
soil_bac_sub_no_bad_Filtered_final_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/soil_bac_sub_no_bad_Filtered_final_non_normalized.rds")
rhizo_bac_sub_no_bad_Filtered_final_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/rhizo_bac_sub_no_bad_Filtered_final_non_normalized.rds")
all_samples_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDs/all_samples_non_normalized.rds")


# ======================================
# Assign Palette colorS to phyla groups
# ======================================
phylum_colors_1 <- c(
  "gray", "mediumaquamarine","plum","salmon", "yellowgreen","dodgerblue")

phylum_colors_2 <- c(
  "gray", "mediumaquamarine", "plum", "gold", "lightsteelblue", 
  "salmon", "mediumpurple", "seagreen", "black", 
  "darkslategray", "indianred", "yellowgreen", "dodgerblue", "darkkhaki")



# ==========================================
# Bacterial composition across sample types
# ==========================================

#For all samples types#

all_ab <- all_samples_non_normalized
sample_data(all_ab)$Growing_Degree_Days <- as.factor(sample_data(all_ab)$Growing_Degree_Days)
all_sample <- all_ab

#Get count of phyla
all_count <- table(phyloseq::tax_table(all_sample)[, "Phylum"])
all_count

#transform to relative abundance
phylumabundance <- all_sample %>%
  tax_glom(taxrank = "Phylum") %>%                        # Set to smallest taxonomic level you are interested in
  transform_sample_counts(function(x) {x/sum(x)} ) %>%   # Transform to rel. abundance
  psmelt()                                               # Melt to long format
phylumabundance

#Make a data frame to include selected variables
all <- phylumabundance %>%
  select(Phylum, Sample_Type, Abundance, Growing_Degree_Days) %>%
  filter(Abundance != 0) %>%
  mutate( Phylum = as.character(Phylum))
all

#Make another data frame to include relative abundance
all_phylum_new <- all %>%
  select(Phylum, Sample_Type, Abundance, Growing_Degree_Days) %>%  #choose variables to work with
  group_by(Sample_Type) %>%                   #group by variables used to plot NOT taxonomic level
  mutate(totalSum = sum(Abundance)) %>%                                 #calculate total abundance of each Phylum
  ungroup() %>%                                                         #remove grouping variables
  group_by(Sample_Type, Phylum) %>%           #group by same variables PLUS taxonomic level
  summarise(                                                            
    Abundance = sum(Abundance),                                         #sum abundance in each phylum for each  group
    totalSum,
    RelAb = Abundance/totalSum) %>%                                     #calculate relative abundance
  unique()
all_phylum_new

# find Phyla whose relative abundance is less than 1%
Others <- all_phylum_new[all_phylum_new$RelAb <= 0.01,]$Phylum

# change their name to "Others"
all_phylum_new[all_phylum_new$Phylum %in% Others,]$Phylum <- 'Others'
length(unique(all_phylum_new$Phylum))

#reorder the phyla to assign the color palette
phyla <- c("Others","Acidobacteriota","Actinobacteriota","Chloroflexi",
           "Planctomycetota","Proteobacteria")
all_phylum_new$Phylum <- factor(all_phylum_new$Phylum, levels = c(phyla))
all_phylum_new

write.csv(all_phylum_new, "C:/Aduragbemi/Manuscript/Review 5/phyla abundance relativity/all_count_relativity.csv")


#re-order the sample type accordingly and covert them to factor
samples <- c("Soil","Rhizosphere","Root","Stalk")
all_phylum_new$Sample_Type <- factor(all_phylum_new$Sample_Type, levels = c(samples))


#Make the final plot
combined_final_1 <- ggplot(all_phylum_new)+
  geom_col(mapping = aes(x = Sample_Type , y = RelAb, fill = Phylum),  
           position = "stack",  show.legend = TRUE) +
  ylab("Relative Abundance") +
  xlab("")+
  scale_fill_manual(values = phylum_colors_1) +
  theme_classic() +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15)),
        axis.title = element_text(family = "serif", size = (15), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,size=15, colour = "black", face = "bold"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black", face = "bold"),
        legend.text = element_text(size = 15),
        legend.title = element_text(face="bold", size = 15)) +
  theme(legend.position="right") +
  guides(fill=guide_legend(title = "Phylum")) +
  ggtitle("")

combined_final_1

#save as image
ggsave("combined_final_1.png", plot = combined_final_1 , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("combined_final_1.pdf", plot = combined_final_1, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")



# ==================================================================
# Bacterial composition across sample types and growing degree days
# ==================================================================

#For soil samples#

#re-assign variables and convert GDD to a factor
soil_ab <- soil_bac_sub_no_bad_Filtered_final_non_normalized
sample_data(soil_ab)$Growing_Degree_Days <- as.factor(sample_data(soil_ab)$Growing_Degree_Days)
soil_sample <- soil_ab

#Get count of phyla
soil_count <- table(phyloseq::tax_table(soil_sample)[, "Phylum"])

#transform to relative abundance
phylumabundance_soil <- soil_sample %>%
  tax_glom(taxrank = "Phylum") %>%                        # Set to smallest taxonomic level of interest
  transform_sample_counts(function(x) {x/sum(x)} ) %>%   # Transform to rel. abundance
  psmelt()                                               # Melt to long format

phylumabundance_soil


#Make a data frame to include selected variables
all_soil <- phylumabundance_soil %>%
  select(Phylum, Sample_Type, Abundance, Growing_Degree_Days) %>%
  filter(Abundance != 0) %>%
  mutate( Phylum = as.character(Phylum))

all_soil

#Make another data frame to include relative abundance
soil_phylum <- all_soil %>%
  select(Phylum, Sample_Type, Abundance, Growing_Degree_Days) %>%  #choose variables to work with
  group_by(Growing_Degree_Days) %>%                   #group by variables used to plot NOT taxonomic level
  mutate(totalSum = sum(Abundance)) %>%                                 #calculate total abundance of each Phylum
  ungroup() %>%                                                         #remove grouping variables
  group_by(Growing_Degree_Days, Phylum, Sample_Type) %>%           #group by same variables PLUS taxonomic level
  summarise(                                                            
    Abundance = sum(Abundance),                                         #sum abundance in each phylum for each  group
    totalSum,
    RelAb = Abundance/totalSum) %>%                                     #calculate relative abundance
  unique()

soil_phylum

# find Phyla whose relative abundunce is less than 1%
Others <- soil_phylum[soil_phylum$RelAb <= 0.01,]$Phylum

# change their name to "Others"
soil_phylum[soil_phylum$Phylum %in% Others,]$Phylum <- 'Others'
length(unique(soil_phylum$Phylum))

#write the soil_phylum to excel sheet
write.csv(soil_phylum, "C:/Aduragbemi/Manuscript/Review 5/phyla abundance relativity/soil_count_relativity.csv")



#For Rhizosphere samples#

#re-assign variables and convert GDD to a factor
rhizo_ab <- rhizo_bac_sub_no_bad_Filtered_final_non_normalized
sample_data(rhizo_ab)$Growing_Degree_Days <- as.factor(sample_data(rhizo_ab)$Growing_Degree_Days)
rhizo_sample <- rhizo_ab

#Get count of phyla
rhizo_count <- table(phyloseq::tax_table(rhizo_sample)[, "Phylum"])

#transform to relative abundance
phylumabundance_rhizo <- rhizo_sample %>%
  tax_glom(taxrank = "Phylum") %>%                        # Set to smallest taxonomic level of interest
  transform_sample_counts(function(x) {x/sum(x)} ) %>%   # Transform to rel. abundance
  psmelt()                                               # Melt to long format

phylumabundance_rhizo


#Make a data frame to include selected variables
all_rhizo <- phylumabundance_rhizo %>%
  select(Phylum, Sample_Type, Abundance, Growing_Degree_Days) %>%
  filter(Abundance != 0) %>%
  mutate( Phylum = as.character(Phylum))

all_rhizo

#Make another data frame to include relative abundance
rhizo_phylum <- all_rhizo %>%
  select(Phylum, Sample_Type, Abundance, Growing_Degree_Days) %>%  #choose variables to work with
  group_by(Growing_Degree_Days) %>%                   #group by variables used to plot NOT taxonomic level
  mutate(totalSum = sum(Abundance)) %>%                                 #calculate total abundance of each Phylum
  ungroup() %>%                                                         #remove grouping variables
  group_by(Growing_Degree_Days, Phylum, Sample_Type) %>%           #group by same variables PLUS taxonomic level
  summarise(                                                            
    Abundance = sum(Abundance),                                         #sum abundance in each phylum for each  group
    totalSum,
    RelAb = Abundance/totalSum) %>%                                     #calculate relative abundance
  unique()

rhizo_phylum

# find Phyla whose relative abundunce is less than 1%
Others <- rhizo_phylum[rhizo_phylum$RelAb <= 0.01,]$Phylum

# change their name to "Others"
rhizo_phylum[rhizo_phylum$Phylum %in% Others,]$Phylum <- 'Others'
length(unique(rhizo_phylum$Phylum))

#write the rhizo_phylum to excel sheet
write.csv(rhizo_phylum, "C:/Aduragbemi/Manuscript/Review 5/phyla abundance relativity/rhizo_count_relativity.csv")


#For Root samples#

#re-assign variables and convert GDD to a factor
root_ab <- root_bac_sub_no_bad_Filtered_final_non_normalized
sample_data(root_ab)$Growing_Degree_Days <- as.factor(sample_data(root_ab)$Growing_Degree_Days)
root_sample <- root_ab

#Get count of phyla
root_count <- table(phyloseq::tax_table(root_sample)[, "Phylum"])

#transform to relative abundance
phylumabundance_root <- root_sample %>%
  tax_glom(taxrank = "Phylum") %>%                        # Set to smallest taxonomic level of interest
  transform_sample_counts(function(x) {x/sum(x)} ) %>%   # Transform to rel. abundance
  psmelt()                                               # Melt to long format

phylumabundance_root


#Make a data frame to include selected variables
all_root <- phylumabundance_root %>%
  select(Phylum, Sample_Type, Abundance, Growing_Degree_Days) %>%
  filter(Abundance != 0) %>%
  mutate( Phylum = as.character(Phylum))

all_root

#Make another data frame to include relative abundance
root_phylum <- all_root %>%
  select(Phylum, Sample_Type, Abundance, Growing_Degree_Days) %>%  #choose variables to work with
  group_by(Growing_Degree_Days) %>%                   #group by variables used to plot NOT taxonomic level
  mutate(totalSum = sum(Abundance)) %>%                                 #calculate total abundance of each Phylum
  ungroup() %>%                                                         #remove grouping variables
  group_by(Growing_Degree_Days, Phylum, Sample_Type) %>%           #group by same variables PLUS taxonomic level
  summarise(                                                            
    Abundance = sum(Abundance),                                         #sum abundance in each phylum for each  group
    totalSum,
    RelAb = Abundance/totalSum) %>%                                     #calculate relative abundance
  unique()

root_phylum

# find Phyla whose relative abundunce is less than 1%
Others <- root_phylum[root_phylum$RelAb <= 0.01,]$Phylum

# change their name to "Others"
root_phylum[root_phylum$Phylum %in% Others,]$Phylum <- 'Others'
length(unique(root_phylum$Phylum))

#write the root_phylum to excel sheet
write.csv(root_phylum, "C:/Aduragbemi/Manuscript/Review 5/phyla abundance relativity/root_count_relativity.csv")



#For Stalk samples#

#re-assign variables and convert GDD to a factor
stalk_ab <- stalk_bac_sub_no_bad_Filtered_final_non_normalized
sample_data(stalk_ab)$Growing_Degree_Days <- as.factor(sample_data(stalk_ab)$Growing_Degree_Days)
stalk_sample <- stalk_ab

#Get count of phyla
stalk_count <- table(phyloseq::tax_table(stalk_sample)[, "Phylum"])

#transform to relative abundance
phylumabundance_stalk <- stalk_sample %>%
  tax_glom(taxrank = "Phylum") %>%                        # Set to smallest taxonomic level of interest
  transform_sample_counts(function(x) {x/sum(x)} ) %>%   # Transform to rel. abundance
  psmelt()                                               # Melt to long format

phylumabundance_stalk


#Make a data frame to include selected variables
all_stalk <- phylumabundance_stalk %>%
  select(Phylum, Sample_Type, Abundance, Growing_Degree_Days) %>%
  filter(Abundance != 0) %>%
  mutate( Phylum = as.character(Phylum))

all_stalk

#Make another data frame to include relative abundance
stalk_phylum <- all_stalk %>%
  select(Phylum, Sample_Type, Abundance, Growing_Degree_Days) %>%  #choose variables to work with
  group_by(Growing_Degree_Days) %>%                   #group by variables used to plot NOT taxonomic level
  mutate(totalSum = sum(Abundance)) %>%                                 #calculate total abundance of each Phylum
  ungroup() %>%                                                         #remove grouping variables
  group_by(Growing_Degree_Days, Phylum, Sample_Type) %>%           #group by same variables PLUS taxonomic level
  summarise(                                                            
    Abundance = sum(Abundance),                                         #sum abundance in each phylum for each  group
    totalSum,
    RelAb = Abundance/totalSum) %>%                                     #calculate relative abundance
  unique()

stalk_phylum

# find Phyla whose relative abundunce is less than 1%
Others <- stalk_phylum[stalk_phylum$RelAb <= 0.01,]$Phylum

# change their name to "Others"
stalk_phylum[stalk_phylum$Phylum %in% Others,]$Phylum <- 'Others'
length(unique(stalk_phylum$Phylum))

#write the stalk_phylum to excel sheet
write.csv(stalk_phylum, "C:/Aduragbemi/Manuscript/Review 5/phyla abundance relativity/stalk_count_relativity.csv")


#Next is to combine all tissue type to form a data frame
combined_samples <- bind_rows(stalk_phylum, rhizo_phylum, root_phylum, soil_phylum)
combined_samples

#reorder the phyla to assign the color palette
phyla <- c("Others","Acidobacteriota","Actinobacteriota","Armatimonadota","Bacteroidota",
           "Chloroflexi","Firmicutes","Gemmatimonadota","Latescibacterota","Methylomirabilota",
           "Myxococcota","Planctomycetota","Proteobacteria","Verrucomicrobiota")

#make phylum variable a factor
combined_samples$Phylum <- factor(combined_samples$Phylum, levels = c(phyla))
combined_samples

#re-order the sample type accordingly and covert them to factor
samples <- c("Soil","Rhizosphere","Root","Stalk")
combined_samples$Sample_Type <- factor(combined_samples$Sample_Type, levels = c(samples))

#make the final plot
combined_final <- ggplot(combined_samples)+
  geom_col(mapping = aes(x = Growing_Degree_Days, y = RelAb, fill = Phylum),  
           position = "stack",  show.legend = TRUE) +
  facet_grid(cols = vars(Sample_Type), )+
  ylab("Relative Abundance") +
  xlab("Growing Degree Days")+
  scale_fill_manual(values = phylum_colors_2) +
  theme_classic() +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (10)),
        axis.title = element_text(family = "serif", size = (15), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,size=15, colour = "black"),
        axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,size=15, colour = "black"),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 15)) +
  theme(legend.position="right") +
  guides(colour = guide_legend(override.aes = list(size = 1.2))) +
  ggtitle("")

combined_final

#save as image
ggsave("combined_final.png", plot = combined_final , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")


#save as pdf
ggsave("combined_final.pdf", plot = combined_final, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")


