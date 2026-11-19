#Load Libraries
library(phyloseq)
library(decontam)
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)


# =========================
# read in the raw data RDS
# =========================


#Metadata#
SAMP.bac <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/SAMP.bac.rds")
SAMP.bac

#OTU table#
OTU.bac <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/OTU.bac.rds")
OTU.bac

#Taxonomy#
TAX.bac <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/taxonomy.rds")
TAX.bac 

#Fasta#
FASTA.bac <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/FASTA.bac.rds")
FASTA.bac

#Phylogenetic tree#
tree <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/tree.rds")
tree

#Create Initial Phyloseq object#
bac.unedited <- phyloseq::phyloseq(OTU.bac, TAX.bac, FASTA.bac, SAMP.bac, tree)
bac.unedited


#change the name to ASVs for easy identification
names(FASTA.bac) <- taxa_names(bac.unedited)
bac.unedited <- merge_phyloseq(bac.unedited, FASTA.bac)
taxa_names(bac.unedited) <- paste0("ASV", seq(ntaxa(bac.unedited)))
bac.unedited@tax_table

#check if the taxa_names of phyloseq and rownames for otu table are the same
identical(taxa_names(bac.unedited),rownames(bac.unedited@otu_table))



# ========================
# Perform some  filtering
# ========================

#Taxonomy filtering#

#remove ASVs that are mitochondria, chloroplast, or Unassigned at the kingdom level 
bac_no_chloro <- bac.unedited %>% 
  phyloseq::subset_taxa(Order != "Chloroplast") %>%
  phyloseq::subset_taxa(Family != "Mitochondria") %>%
  phyloseq::subset_taxa(Kingdom == "Bacteria") %>%
  phyloseq::subset_taxa(Phylum != "Unclassified ") %>%
  phyloseq::subset_samples(Sample_Type != "Inrow_Soil")


#next remove sample with less than 500 reads
bac_no_chloro <- prune_samples(sample_sums(bac_no_chloro) > 500, bac_no_chloro)
bac_no_chloro

# Check total taxa in each sample again
sample_sums(bac_no_chloro)
sort(sample_sums(bac_no_chloro))


#Data filtering#

#Convert Growing_Degree_Days to factor 
sample_data(bac_no_chloro)$Growing_Degree_Days <- as.factor(sample_data(bac_no_chloro)$Growing_Degree_Days)
sample_data(bac_no_chloro)$Location <- as.factor(sample_data(bac_no_chloro)$Location)


#fix location names
sample_data(bac_no_chloro)$Location <- gsub("_", " ", sample_data(bac_no_chloro)$Location)
sample_data(bac_no_chloro)$Location <- gsub("Watkinsville IronHorse", "Watkinsville", sample_data(bac_no_chloro)$Location)
sample_data(bac_no_chloro)$Location <- gsub("Plains Production corn", "Plains-Production", sample_data(bac_no_chloro)$Location)
sample_data(bac_no_chloro)$Location <- gsub("Plains SWVT", "Plains-SWVT", sample_data(bac_no_chloro)$Location)
sample_data(bac_no_chloro)$Location <- factor(sample_data(bac_no_chloro)$Location, levels = c("Pavo", "Arlington", "Tifton-Irrigated", "Tifton-Dryland", "Dawson", "Fitzgerald", "Plains-Production", "Plains-SWVT", "Hawkinsville", "Fort Valley", "Midville", "Tennile", "Wadley", "Athens", "Watkinsville", "Cave Spring", "Rome", "Blairsville"))
sample_data(bac_no_chloro)$Location
sample_data(bac_no_chloro)$RepName <- gsub("[ab]", "", sample_data(bac_no_chloro)$SampleID, perl=TRUE)

bac_no_chloro_all <- bac_no_chloro 


#Filtering the whole dataset/samples
bac_no_chloro_all <- prune_taxa(taxa_sums(bac_no_chloro_all)>1, bac_no_chloro_all)
bac_no_chloro_all



# =========================================
# Removing contaminant / decontaminate#
# =========================================

sample_data(bac_no_chloro_all)$is.neg <- sample_data(bac_no_chloro_all)$Sample_or_Control == "Control"
all_contamdf.prev <- isContaminant(bac_no_chloro_all, method="prevalence", neg="is.neg", threshold = 0.5)
all_badTaxa <- rownames(all_contamdf.prev[all_contamdf.prev$contaminant == TRUE,])

print(all_badTaxa) # 717

#calculate the Proportion of ASVs identified as contaminants for each sample type
length(which(all_contamdf.prev$contaminant)) / nrow(all_contamdf.prev)   # 0.02092148

#continue with decontaminate
all_ps.pa <- transform_sample_counts(bac_no_chloro_all, function(abund) 1*(abund>0))
all_ps.pa.neg <- prune_samples(sample_data(all_ps.pa)$Sample_or_Control == "Control", all_ps.pa)
all_ps.pa.pos <- prune_samples(sample_data(all_ps.pa)$Sample_or_Control == "Sample", all_ps.pa)

# Make data frame of prevalence in positive and negative samples
all_df.pa <- data.frame(pa.pos=taxa_sums(all_ps.pa.pos), pa.neg=taxa_sums(all_ps.pa.neg),
                          contaminant=all_contamdf.prev$contaminant)
all_df.pa

#make a plot
all_decontaminate.bac <- ggplot(data=all_df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
  geom_point(size =10) +
  xlab("Prevalence (Negative Controls)") + 
  ylab("Prevalence (True Samples)") +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (40), colour = "black"),
        axis.title = element_text(family = "serif", size = (40), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=0,hjust=0.5,vjust=1.5, size=30, face = "bold", colour = "black"),
        axis.text.y=element_text(angle=0,hjust=1.1,vjust=0.4, size = 30,  face = "bold", colour = "black"),
        legend.title = element_text(colour = "black", size = 30, face = "bold"),
        legend.text = element_text(size = 30, colour = "black"))  +
  ggtitle("All Sample Types")

all_decontaminate.bac

#save the figure
#save as image
ggsave("all_decontaminate.bac.png", plot = all_decontaminate.bac, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("all_decontaminate.bac.pdf", plot = all_decontaminate.bac, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 10, height = 5, units = c("in"), device = "pdf")


#remove the contaminants and save the data free from contaminants as
all_goodTaxa <- setdiff(taxa_names(bac_no_chloro_all), all_badTaxa)
all_samples_no_bad_taxa <- prune_taxa(all_goodTaxa, bac_no_chloro_all)
all_samples_no_bad_taxa 
sort(sample_sums(all_samples_no_bad_taxa))

all_samples_non_normalized <- all_samples_no_bad_taxa  
all_samples_non_normalized
sort(sample_sums(all_samples_non_normalized))

#remove the blanks before saving to RDS
all_samples_non_normalized <- subset_samples(all_samples_non_normalized, Sample_Type == "Rhizosphere" | Sample_Type == "Root" | Sample_Type == "Soil" | Sample_Type == "Stalk")
all_samples_non_normalized

#Save in the non normalized ASVs for all sample type together here
saveRDS(all_samples_non_normalized, file = "C:/Aduragbemi/Microbiome/RDs/all_samples_non_normalized.rds")

#Read in the non normalized ASVs for all sample type together here
all_samples_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDs/all_samples_non_normalized.rds")
all_samples_non_normalized


# Check total taxa in each sample again
sample_sums(all_samples_non_normalized)
sort(sample_sums(all_samples_non_normalized))


#subset the samples to only tissue removing the blanks#
all_samples_non_normalized_final <- subset_samples(all_samples_non_normalized, Sample_Type == "Root" | Sample_Type == "Soil" | Sample_Type == "Stalk" | Sample_Type =="Rhizosphere")
all_samples_non_normalized_final


#now save the normalized OTUs for all sample type together here
saveRDS(all_samples_non_normalized_final, file = "C:/Aduragbemi/Microbiome/RDS/all_samples_non_normalized_final.rds")



#Subset the decontaminated phyloseq "all_samples_non_normalized"#

#subset data by sample type
bac_no_chloro_rhizo <- subset_samples(all_samples_non_normalized, Sample_Type == "Rhizosphere" | Sample_Type == "blank-Rhizosphere")
bac_no_chloro_rhizo
bac_no_chloro_root <- subset_samples(all_samples_non_normalized, Sample_Type == "Root" | Sample_Type == "blank-Root")
bac_no_chloro_root
bac_no_chloro_stalk <- subset_samples(all_samples_non_normalized, Sample_Type == "Stalk" | Sample_Type == "blank-Stalk")
bac_no_chloro_stalk
bac_no_chloro_soil <- subset_samples(all_samples_non_normalized, Sample_Type == "Soil"| Sample_Type == "blank-Soil" ) 
bac_no_chloro_soil


#prune subset samples
bac_no_chloro_soil <- prune_taxa(taxa_sums(bac_no_chloro_soil)>1, bac_no_chloro_soil)
bac_no_chloro_root <- prune_taxa(taxa_sums(bac_no_chloro_root)>1, bac_no_chloro_root)
bac_no_chloro_rhizo <- prune_taxa(taxa_sums(bac_no_chloro_rhizo)>1, bac_no_chloro_rhizo)
bac_no_chloro_stalk <- prune_taxa(taxa_sums(bac_no_chloro_stalk)>1, bac_no_chloro_stalk)


#subset samples again to include just the growing degree days of 600 and 1400
bac_no_chloro_rhizo <- subset_samples(bac_no_chloro_rhizo , Growing_Degree_Days == "600" | Growing_Degree_Days == "1400" | Growing_Degree_Days == "NA-Rhizosphere")
bac_no_chloro_rhizo
bac_no_chloro_root <- subset_samples(bac_no_chloro_root, Growing_Degree_Days == "600" | Growing_Degree_Days == "1400" | Growing_Degree_Days == "NA-Root")
bac_no_chloro_root
bac_no_chloro_stalk <- subset_samples(bac_no_chloro_stalk, Growing_Degree_Days == "600" | Growing_Degree_Days == "1400" | Growing_Degree_Days == "NA-Stalk")
bac_no_chloro_stalk
bac_no_chloro_soil <- subset_samples(bac_no_chloro_soil, Growing_Degree_Days == "600"| Growing_Degree_Days == "1400" | Growing_Degree_Days == "NA-Soil") 
bac_no_chloro_soil


# ===============
# Data filtering
# ===============


#remove samples  blank that are outliers
bac_no_chloro_rhizo <- subset_samples(bac_no_chloro_rhizo)
bac_no_chloro_rhizo

bac_no_chloro_root <- subset_samples(bac_no_chloro_root, SampleID != "P3BR")
bac_no_chloro_root

bac_no_chloro_soil <- subset_samples(bac_no_chloro_soil, SampleID != "P5BS") 
bac_no_chloro_soil

#Remove samples <500 reads (redone Just to be sure)

rhizo_bac_sub_no_bad_Filtered <- prune_samples(sample_sums(bac_no_chloro_rhizo) > 500, bac_no_chloro_rhizo)
root_bac_sub_no_bad_Filtered <- prune_samples(sample_sums(bac_no_chloro_root) > 500, bac_no_chloro_root)
soil_bac_sub_no_bad_Filtered <- prune_samples(sample_sums(bac_no_chloro_soil) > 500, bac_no_chloro_soil)
stalk_bac_sub_no_bad_Filtered <- prune_samples(sample_sums(bac_no_chloro_stalk) > 500, bac_no_chloro_stalk)


# ================================================================================================================================
# Uncover the variability in number of reads per sample for each sample_type data set
# Make a data frame of sample sums and sample types, then make violin plot to show distribution of sample sums for each sample type
# ================================================================================================================================

#for Rhizosphere samples

rhizo_bac_sub_no_bad_sequencing_read_depth <- data.frame(SampleID = sample_names(rhizo_bac_sub_no_bad_Filtered), 
                                                         rhizo_Depth = sample_sums(rhizo_bac_sub_no_bad_Filtered), 
                                                         rhizo_Tissue = sample_data(rhizo_bac_sub_no_bad_Filtered)$Sample_Type, 
                                                         rhizo_Sample_or_Control = sample_data(rhizo_bac_sub_no_bad_Filtered)$Sample_or_Control)
rhizo_bac_sub_no_bad_sequencing_read_depth

#for Stalk samples
stalk_bac_sub_no_bad_sequencing_read_depth <- data.frame(SampleID = sample_names(stalk_bac_sub_no_bad_Filtered), 
                                                         stalk_Depth = sample_sums(stalk_bac_sub_no_bad_Filtered), 
                                                         stalk_Tissue = sample_data(stalk_bac_sub_no_bad_Filtered)$Sample_Type, 
                                                         stalk_Sample_or_Control = sample_data(stalk_bac_sub_no_bad_Filtered)$Sample_or_Control)
stalk_bac_sub_no_bad_sequencing_read_depth

#for Root samples
root_bac_sub_no_bad_sequencing_read_depth <- data.frame(SampleID = sample_names(root_bac_sub_no_bad_Filtered), 
                                                        root_Depth = sample_sums(root_bac_sub_no_bad_Filtered), 
                                                        root_Tissue = sample_data(root_bac_sub_no_bad_Filtered)$Sample_Type, 
                                                        root_Sample_or_Control = sample_data(root_bac_sub_no_bad_Filtered)$Sample_or_Control)

root_bac_sub_no_bad_sequencing_read_depth

#for Soil samples
soil_bac_sub_no_bad_sequencing_read_depth <- data.frame(SampleID = sample_names(soil_bac_sub_no_bad_Filtered), 
                                                        soil_Depth = sample_sums(soil_bac_sub_no_bad_Filtered), 
                                                        soil_Tissue = sample_data(soil_bac_sub_no_bad_Filtered)$Sample_Type, 
                                                        soil_Sample_or_Control = sample_data(soil_bac_sub_no_bad_Filtered)$Sample_or_Control)
soil_bac_sub_no_bad_sequencing_read_depth


#write out in csv each file edit then and import them back to plot sequencing plot for each tissue#

write.csv(soil_bac_sub_no_bad_sequencing_read_depth, "C:/Aduragbemi/Microbiome/Result files/soil_bac_sub_no_bad_sequencing_read_depth.csv")
write.csv(root_bac_sub_no_bad_sequencing_read_depth,"C:/Aduragbemi/Microbiome/Result files/root_bac_sub_no_bad_sequencing_read_depth.csv")
write.csv(stalk_bac_sub_no_bad_sequencing_read_depth, "C:/Aduragbemi/Microbiome/Result files/stalk_bac_sub_no_bad_sequencing_read_depth.csv")
write.csv(rhizo_bac_sub_no_bad_sequencing_read_depth, "C:/Aduragbemi/Microbiome/Result files/rhizo_bac_sub_no_bad_sequencing_read_depth.csv")


#read all back as a single file all_bac_sub_no_bad_Filtered_sequencing_read_depth
all_samples <- read.csv("C:/Aduragbemi/Microbiome/Result files/all_bac_sub_no_bad_sequencing_read_depth.csv")
all_samples

samples <- c("Soil","Rhizosphere","Root","Stalk")


all_samples$Sample_Type <- factor(all_samples$Sample_Type, levels = c(samples))


#plot the graph for sequencing depth of each sample type

all_samples_plot <-  ggplot(all_samples, 
                            aes(x=Sample_Type, y=Sequencing_Depth)) +
  theme(legend.position = "") +
  geom_violin(width =0.5, aes(fill = Sample_Type), position = position_dodge(1), scale = "width") +
  scale_fill_brewer(palette="Dark2")+
  geom_boxplot(width=0.1, fill= "white")+
  ylab ("Sequencing Depth") +
  xlab ("") +
  theme(plot.title = element_text(family = "serif", face = "bold", size = (15), colour = "black"),
        axis.title = element_text(family = "serif", size = (15), colour = "black", face = "bold"),
        axis.text.x=element_text(angle=0,hjust=0.5,vjust=1.0,size=15, face = "bold", family = "serif", colour = "black"),
        axis.text.y=element_text(angle=0,hjust=1.1,vjust=0.4,size=15, face = "bold", family = "serif", colour = "black"),
        legend.title = element_text(colour = "black", size = 7, face = "bold"),
        legend.text = element_text(size = 7, colour = "black")) +
  ggtitle("Sequencing Depth for all Sample Types")

all_samples_plot

#save the figure
#save as image
ggsave("all_samples_plot.png", plot = all_samples_plot, path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 5, height = 5, units = c("in"), device = "png")

#save as pdf
ggsave("all_samples_plot.pdf", plot = all_samples_plot, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 5, height = 5, units = c("in"), device = "pdf")



###Calculate the reads per sample###

#Soil sample
sample.sums_soil <- data.frame(sample_sums(soil_bac_sub_no_bad_Filtered))

soil_read.dist.bac <- ggplot(sample.sums_soil, aes(x = sample_sums.soil_bac_sub_no_bad_Filtered.)) +
  geom_histogram(color = "black", fill = cbbPalette[[4]]) + 
  theme_classic() +
  xlab("Read Depth") + 
  ggtitle("Soil")
soil_read.dist.bac 

sum(sample_sums(soil_bac_sub_no_bad_Filtered)) # total reads = 1009206
median(sample_sums(soil_bac_sub_no_bad_Filtered)) # 10664.5


#Rhizosphere sample
sample.sums_rhizo <- data.frame(sample_sums(rhizo_bac_sub_no_bad_Filtered))

rhizo_read.dist.bac <- ggplot(sample.sums_rhizo, aes(x = sample_sums.rhizo_bac_sub_no_bad_Filtered.)) +
  geom_histogram(color = "black", fill = cbbPalette[[4]]) + 
  theme_classic() +
  xlab("Read Depth") + 
  ggtitle("Rhizosphere")
rhizo_read.dist.bac 

sum(sample_sums(rhizo_bac_sub_no_bad_Filtered)) # total reads = 1547065
median(sample_sums(rhizo_bac_sub_no_bad_Filtered)) # 12978



#Stalk sample
sample.sums_stalk <- data.frame(sample_sums(stalk_bac_sub_no_bad_Filtered))

stalk_read.dist.bac <- ggplot(sample.sums_stalk, aes(x = sample_sums.stalk_bac_sub_no_bad_Filtered.)) +
  geom_histogram(color = "black", fill = cbbPalette[[4]]) + 
  theme_classic() +
  xlab("Read Depth") + 
  ggtitle("Stalk")
stalk_read.dist.bac 

sum(sample_sums(stalk_bac_sub_no_bad_Filtered)) # total reads = 310016
median(sample_sums(stalk_bac_sub_no_bad_Filtered)) # 2602.5



#Root sample
sample.sums_root <- data.frame(sample_sums(root_bac_sub_no_bad_Filtered))

root_read.dist.bac <- ggplot(sample.sums_root, aes(x = sample_sums.root_bac_sub_no_bad_Filtered.)) +
  geom_histogram(color = "black", fill = cbbPalette[[4]]) + 
  theme_classic() +
  xlab("Read Depth") + 
  ggtitle("Root")
root_read.dist.bac 

sum(sample_sums(root_bac_sub_no_bad_Filtered)) # total reads = 788170
median(sample_sums(root_bac_sub_no_bad_Filtered)) # 4881


#plot the  graph

all_read.dist.bac <- ggarrange(soil_read.dist.bac, stalk_read.dist.bac, root_read.dist.bac, rhizo_read.dist.bac, ncol = 2, nrow = 2, labels = c("A","B","C","D"))
all_read.dist.bac # Figure for finaplot

#save the figure
#save as image
ggsave("all_read.dist.bac.png", plot = all_read.dist.bac, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 30, height = 15, units = c("in"), device = "png")

#save as pdf
ggsave("all_read.dist.bac.pdf", plot = all_read.dist.bac, path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 30, height = 15, units = c("in"), device = "pdf")



#final filtering by subsetting samples before saving it to RDS#
soil_bac_sub_no_bad_Filtered_final_non_normalized <- subset_samples(soil_bac_sub_no_bad_Filtered, Sample_Type == "Soil") 
root_bac_sub_no_bad_Filtered_final_non_normalized <- subset_samples(root_bac_sub_no_bad_Filtered, Sample_Type == "Root")
rhizo_bac_sub_no_bad_Filtered_final_non_normalized <- subset_samples(rhizo_bac_sub_no_bad_Filtered, Sample_Type == "Rhizosphere")
stalk_bac_sub_no_bad_Filtered_final_non_normalized <- subset_samples(stalk_bac_sub_no_bad_Filtered, Sample_Type == "Stalk")



#Save as RDS
saveRDS(root_bac_sub_no_bad_Filtered_final_non_normalized, file = "C:/Aduragbemi/Microbiome/RDS/root_bac_sub_no_bad_Filtered_final_non_normalized.rds")
saveRDS(stalk_bac_sub_no_bad_Filtered_final_non_normalized, file = "C:/Aduragbemi/Microbiome/RDS/stalk_bac_sub_no_bad_Filtered_final_non_normalized.rds")
saveRDS(soil_bac_sub_no_bad_Filtered_final_non_normalized, file = "C:/Aduragbemi/Microbiome/RDS/soil_bac_sub_no_bad_Filtered_final_non_normalized.rds")
saveRDS(rhizo_bac_sub_no_bad_Filtered_final_non_normalized, file = "C:/Aduragbemi/Microbiome/RDS/rhizo_bac_sub_no_bad_Filtered_final_non_normalized.rds")
