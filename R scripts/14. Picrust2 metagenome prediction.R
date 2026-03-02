#Load Libraries
library(Biostrings)
library(magrittr)
library(phyloseq)
library(biomformat)
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(ggplot2)
library(ggplotify)
library("KEGGREST")


# ===============================
# Load in color for sample types
# ===============================

sample_colors <- c(
  "Soil"        = "#008080",   
  "Rhizosphere" = "#D2691E",   
  "Root"        = "#6A5ACD",   
  "Stalk"       = "#FF1493"    
)

# ===============================
# read in the non-normalized data
# ===============================

all_picrust <- readRDS(file = "C:/Aduragbemi/Microbiome/RDs/all_samples_non_normalized.rds")
root_picrust <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/root_bac_sub_no_bad_Filtered_final_non_normalized.rds")
stalk_picrust <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/stalk_bac_sub_no_bad_Filtered_final_non_normalized.rds")
soil_picrust <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/soil_bac_sub_no_bad_Filtered_final_non_normalized.rds")
rhizo_picrust <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/rhizo_bac_sub_no_bad_Filtered_final_non_normalized.rds")


#since i would like to compare the sample types within each GDD, I will subdivide all_picrust into two

#subset samples again to include just the growing degree days of 600 and 1400
all_picrust_600 <- subset_samples(all_picrust , Growing_Degree_Days == "600")
all_picrust_600

all_picrust_1400 <- subset_samples(all_picrust , Growing_Degree_Days == "1400")
all_picrust_1400


# ===============================================
# first extract the refseq for each samples types
# ===============================================

#for root sample

root_picrust %>%
  refseq() %>%
  Biostrings::writeXStringSet("C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/root_asv.fna", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")

#for stalk sample
stalk_picrust %>%
  refseq() %>%
  Biostrings::writeXStringSet("C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/stalk_asv.fna", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")

#for rhizo sample
rhizo_picrust %>%
  refseq() %>%
  Biostrings::writeXStringSet("C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/rhizo_asv.fna", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")

#for soil sample
soil_picrust %>%
  refseq() %>%
  Biostrings::writeXStringSet("C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/soil_asv.fna", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")

#all_picrust_600
all_picrust_600 %>%
  refseq() %>%
  Biostrings::writeXStringSet("C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/all_picrust_600.fna", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")

#all_picrust_1400
all_picrust_1400 %>%
  refseq() %>%
  Biostrings::writeXStringSet("C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/all_picrust_1400.fna", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")



# ================================================================
# Next is to extract the ASVs of each sample types into biom files
# ================================================================

#For soil samples
#check if taxa are rows
taxa_are_rows(soil_picrust)

#if taxa_are_rows=TRUE
soil_otu<-as(otu_table(soil_picrust),"matrix")
soil_otu_biom<-make_biom(data=soil_otu)
write_biom(soil_otu_biom,"C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/soil_otu_biom.biom")

#export otu table as a text file
write.table(otu_table(soil_picrust), "C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/soil_otu.txt",sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)


#for rhizosphere samples

#check if taxa are rows
taxa_are_rows(rhizo_picrust)

#if taxa_are_rows=TRUE
rhizo_otu<-as(otu_table(rhizo_picrust),"matrix")
rhizo_otu_biom<-make_biom(data=rhizo_otu)
write_biom(rhizo_otu_biom,"C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/rhizo_otu_biom.biom")

#export otu table as a text file
write.table(otu_table(rhizo_picrust), "C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/rhizo_otu.txt",sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)


#for root samples

#check if taxa are rows
taxa_are_rows(root_picrust)

#if taxa_are_rows=TRUE
root_otu<-as(otu_table(root_picrust),"matrix")
root_otu_biom<-make_biom(data=root_otu)
write_biom(root_otu_biom,"C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/root_otu_biom.biom")

#export otu table as a text file
write.table(otu_table(root_picrust), "C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/root_otu.txt",sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)


#For stalk samples

#check if taxa are rows
taxa_are_rows(stalk_picrust)

#if taxa_are_rows=TRUE
stalk_otu<-as(otu_table(stalk_picrust),"matrix")
stalk_otu_biom<-make_biom(data=stalk_otu)
write_biom(stalk_otu_biom,"C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/stalk_otu_biom.biom")

#export otu table as a text file
write.table(otu_table(stalk_picrust), "C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/stalk_otu.txt",sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)


#all_picrust_600

#check if taxa are rows
taxa_are_rows(all_picrust_600)

#if taxa_are_rows=TRUE
all_picrust_600_otu<-as(otu_table(all_picrust_600),"matrix")
all_picrust_600_otu_biom<-make_biom(data=all_picrust_600_otu)
write_biom(all_picrust_600_otu_biom,"C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/all_picrust_600_otu_biom.biom")

#export otu table as a text file
write.table(otu_table(all_picrust_600), "C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/all_picrust_600_otu.txt",sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)


#all_picrust_1400

#check if taxa are rows
taxa_are_rows(all_picrust_1400)

#if taxa_are_rows=TRUE
all_picrust_1400_otu<-as(otu_table(all_picrust_1400),"matrix")
all_picrust_1400_otu_biom<-make_biom(data=all_picrust_1400_otu)
write_biom(all_picrust_1400_otu_biom,"C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/all_picrust_1400_otu_biom.biom")

#export otu table as a text file
write.table(otu_table(all_picrust_1400), "C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/all_picrust_1400_otu.txt",sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)


#Next is to extract the metadata for each sample types
#For soil samples
soil_metadata_df <- data.frame(
  slot(soil_picrust, "sam_data"),
  check.names = FALSE
)

write.csv(soil_metadata_df, "C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/soil_metadata.csv")


#For Rhizosphere samples
rhizo_metadata_df <- data.frame(
  slot(rhizo_picrust, "sam_data"),
  check.names = FALSE
)

write.csv(rhizo_metadata_df, "C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/rhizo_metadata.csv")


#For root samples
root_metadata_df <- data.frame(
  slot(root_picrust, "sam_data"),
  check.names = FALSE
)

write.csv(root_metadata_df, "C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/root_metadata.csv")


#For Stalk samples
stalk_metadata_df <- data.frame(
  slot(stalk_picrust, "sam_data"),
  check.names = FALSE
)

write.csv(stalk_metadata_df, "C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/stalk_metadata.csv")


#For all_picrust_600
all_picrust_600_metadata_df <- data.frame(
  slot(all_picrust_600, "sam_data"),
  check.names = FALSE
)

write.csv(all_picrust_600_metadata_df, "C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/all_picrust_600_metadata.csv")



#For all_picrust_1400
all_picrust_1400_metadata_df <- data.frame(
  slot(all_picrust_1400, "sam_data"),
  check.names = FALSE
)

write.csv(all_picrust_1400_metadata_df, "C:/Aduragbemi/Manuscript/Review 5/picrust/raw files/all_picrust_1400_metadata.csv")





# ======================================================================
# PICRUSt2 Functional Analysis Across the Two GDD for each sample types
# ======================================================================

#For Soil sample

#no statistically significant features

# For Rhizosphere sample

#read metadata
metadata_rhizo <- read.csv("C:/Aduragbemi/Manuscript/Review 5/picrust/input/rhizo/rhizo_metadata.csv")

# Make Growing_Degree_Days a factor
metadata_rhizo$Growing_Degree_Days <- as.factor(metadata_rhizo$Growing_Degree_Days)

# Verify the grouping variable
print(table(metadata_rhizo$Growing_Degree_Days))

# Convert KO abundance to KEGG pathway abundance
kegg_abundance_rhizo <- ko2kegg_abundance("C:/Aduragbemi/Manuscript/Review 5/picrust/input/rhizo/final_rhizo.tsv") 

# Perform pathway DAA using DESeq2
daa_results_df_rhizo <- pathway_daa(
    abundance = kegg_abundance_rhizo, 
    metadata = metadata_rhizo, 
    group = "Growing_Degree_Days", 
    daa_method = "DESeq2", select = NULL, 
    reference = NULL)

# Filter for DESeq2 results
daa_sub_method_results_df_rhizo <- daa_results_df_rhizo[daa_results_df_rhizo$method == "DESeq2", ]

# Annotate pathway results
daa_annotated_sub_method_results_df_rhizo <- pathway_annotation(
    pathway = "KO", daa_results_df = 
    daa_sub_method_results_df_rhizo, 
    ko_to_kegg = TRUE
)

# Remove ko ids related to some pathways not implicated in maize microbiome
daa_human_euk_rhizo <- daa_annotated_sub_method_results_df_rhizo %>%
  dplyr::filter(grepl("^Human Diseases|^Organismal Systems|^Genetic Information Processing|^Cellular Processes", pathway_class))

# Extract KO IDs
ko_human_euk_rhizo <- sort(unique(daa_human_euk_rhizo$feature))
ko_human_euk_rhizo

# Remove KO IDs
daa_annotated_sub_method_results_df_rhizo <- daa_annotated_sub_method_results_df_rhizo %>%
  filter(!feature %in% ko_human_euk_rhizo)


# Generate pathway error bar plot
p_rhizo_picrust <- pathway_errorbar(
  abundance = kegg_abundance_rhizo,
  daa_results_df = daa_annotated_sub_method_results_df_rhizo,
  Group = metadata_rhizo$Growing_Degree_Days,
  p_values_threshold = 0.05,
  order = "pathway_class",
  select = NULL,
  ko_to_kegg = TRUE,
  p_value_bar = TRUE,
  colors = NULL,
  x_lab = "pathway_name"
)

p_rhizo_picrust

#save as image
ggsave("p_rhizo_picrust.png", plot = p_rhizo_picrust , 
       path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 45, height = 15, units = c("in"), device = "png")

#save as pdf
ggsave("p_rhizo_picrust.pdf", plot = p_rhizo_picrust, 
       path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 45, height = 15, units = c("in"), device = "pdf")


# Save DAA results
write_csv(
  daa_annotated_sub_method_results_df_rhizo,
  file.path("C:/Aduragbemi/Manuscript/Review 5/picrust/output/daa_annotated_results_DESeq2_rhizo.csv")
)



#For Root sample

#read metadata
metadata_root <- read.csv("C:/Aduragbemi/Manuscript/Review 5/picrust/input/root/root_metadata.csv")

# Make Growing_Degree_Days a factor
metadata_root$Growing_Degree_Days <- as.factor(metadata_root$Growing_Degree_Days)

# Verify the grouping variable
print(table(metadata_root$Growing_Degree_Days))

# Convert KO abundance to KEGG pathway abundance
kegg_abundance_root <- ko2kegg_abundance("C:/Aduragbemi/Manuscript/Review 5/picrust/input/root/final_root.tsv") 

# Perform pathway DAA using DESeq2
daa_results_df_root <- pathway_daa(
  abundance = kegg_abundance_root, 
  metadata = metadata_root, 
  group = "Growing_Degree_Days", 
  daa_method = "DESeq2", select = NULL, 
  reference = NULL)

# Filter for DESeq2 results
daa_sub_method_results_df_root <- daa_results_df_root[daa_results_df_root$method == "DESeq2", ]

# Annotate pathway results
daa_annotated_sub_method_results_df_root <- pathway_annotation(
  pathway = "KO", daa_results_df = 
    daa_sub_method_results_df_root, 
  ko_to_kegg = TRUE
)

# Remove ko ids related to some pathways not implicated in maize microbiome
daa_human_euk_root <- daa_annotated_sub_method_results_df_root %>%
  dplyr::filter(grepl("^Human Diseases|^Organismal Systems|^Genetic Information Processing|^Cellular Processes", pathway_class))

# Extract KO IDs
ko_human_euk_root <- sort(unique(daa_human_euk_root$feature))
ko_human_euk_root

# Remove KO IDs
daa_annotated_sub_method_results_df_root <- daa_annotated_sub_method_results_df_root %>%
  filter(!feature %in% ko_human_euk_root)

# Generate pathway error bar plot
p_root_picrust <- pathway_errorbar(
  abundance = kegg_abundance_root,
  daa_results_df = daa_annotated_sub_method_results_df_root,
  Group = metadata_root$Growing_Degree_Days,
  p_values_threshold = 0.05,
  order = "pathway_class",
  select = NULL,
  ko_to_kegg = TRUE,
  p_value_bar = TRUE,
  colors = NULL,
  x_lab = "pathway_name"
)

p_root_picrust

#save as image
ggsave("p_root_picrust.png", plot = p_root_picrust , 
       path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 45, height = 15, units = c("in"), device = "png")

#save as pdf
ggsave("p_root_picrust.pdf", plot = p_root_picrust, 
       path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 45, height = 15, units = c("in"), device = "pdf")


# Save DAA results
write_csv(
  daa_annotated_sub_method_results_df_root,
  file.path("C:/Aduragbemi/Manuscript/Review 5/picrust/output/daa_annotated_results_DESeq2_root.csv")
)



#For Stalk sample

#read metadata
metadata_stalk <- read.csv("C:/Aduragbemi/Manuscript/Review 5/picrust/input/stalk/stalk_metadata.csv")

# Make Growing_Degree_Days a factor
metadata_stalk$Growing_Degree_Days <- as.factor(metadata_stalk$Growing_Degree_Days)

# Verify the grouping variable
print(table(metadata_stalk$Growing_Degree_Days))

# Convert KO abundance to KEGG pathway abundance
kegg_abundance_stalk <- ko2kegg_abundance("C:/Aduragbemi/Manuscript/Review 5/picrust/input/stalk/final_stalk.tsv") 

# Perform pathway DAA using DESeq2
daa_results_df_stalk <- pathway_daa(
  abundance = kegg_abundance_stalk, 
  metadata = metadata_stalk, 
  group = "Growing_Degree_Days", 
  daa_method = "DESeq2", select = NULL, 
  reference = NULL)

# Filter for DESeq2 results
daa_sub_method_results_df_stalk <- daa_results_df_stalk[daa_results_df_stalk$method == "DESeq2", ]

# Annotate pathway results
daa_annotated_sub_method_results_df_stalk <- pathway_annotation(
  pathway = "KO", daa_results_df = 
    daa_sub_method_results_df_stalk, 
  ko_to_kegg = TRUE,
)

# Remove ko ids related to some pathways not implicated in maize microbiome
daa_human_euk_stalk <- daa_annotated_sub_method_results_df_stalk %>%
  dplyr::filter(grepl("^Human Diseases|^Organismal Systems|^Genetic Information Processing|^Cellular Processes", pathway_class))

# Extract KO IDs
ko_human_euk_stalk <- sort(unique(daa_human_euk_stalk$feature))
ko_human_euk_stalk

# Remove KO IDs
daa_annotated_sub_method_results_df_stalk <- daa_annotated_sub_method_results_df_stalk %>%
  filter(!feature %in% ko_human_euk_stalk)

# Generate pathway error bar plot
p_stalk_picrust <- pathway_errorbar(
  abundance = kegg_abundance_stalk,
  daa_results_df = daa_annotated_sub_method_results_df_stalk,
  Group = metadata_stalk$Growing_Degree_Days,
  p_values_threshold = 0.05,
  order = "pathway_class",
  select = NULL,
  ko_to_kegg = TRUE,
  p_value_bar = TRUE,
  colors = NULL,
  x_lab = "pathway_name"
)

p_stalk_picrust

#save as image
ggsave("p_stalk_picrust.png", plot = p_stalk_picrust , 
       path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 45, height = 15, units = c("in"), device = "png")

#save as pdf
ggsave("p_stalk_picrust.pdf", plot = p_stalk_picrust, 
       path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 45, height = 15, units = c("in"), device = "pdf")


# Save DAA results
write_csv(
  daa_annotated_sub_method_results_df_stalk,
  file.path("C:/Aduragbemi/Manuscript/Review 5/picrust/output/daa_annotated_results_DESeq2_stalk.csv")
)





# =====================================================
# PICRUSt2 Functional Analysis Across the sample types
# =====================================================

#For GDD_600

#read metadata
metadata_picrust_600 <- read.csv("C:/Aduragbemi/Manuscript/Review 5/picrust/input/picrust_600/all_picrust_600_metadata.csv")

# Make Growing_Degree_Days a factor
metadata_picrust_600$Sample_Type <- as.factor(metadata_picrust_600$Sample_Type)

# Verify the grouping variable
print(table(metadata_picrust_600$Sample_Type))

# Convert KO abundance to KEGG pathway abundance
kegg_abundance_picrust_600 <- ko2kegg_abundance("C:/Aduragbemi/Manuscript/Review 5/picrust/input/picrust_600/all_picrust_600.tsv") 

# Perform pathway DAA using DESeq2
daa_results_df_picrust_600 <- pathway_daa(
  abundance = kegg_abundance_picrust_600, 
  metadata = metadata_picrust_600, 
  group = "Sample_Type", 
  daa_method = "DESeq2", select = NULL, 
  reference = NULL)

# Filter for DESeq2 results
daa_sub_method_results_df_picrust_600 <- daa_results_df_picrust_600[daa_results_df_picrust_600$method == "DESeq2", ]

# Annotate pathway results
daa_annotated_sub_method_results_df_picrust_600 <- pathway_annotation(
  pathway = "KO", daa_results_df = 
    daa_sub_method_results_df_picrust_600, 
  ko_to_kegg = TRUE
)

#Keep only significant KO IDs (adjusted p > 1e-3) ie highly significant
daa_600_filtered <- daa_annotated_sub_method_results_df_picrust_600 %>%
  dplyr::filter(
    !grepl("^Human Diseases|^Organismal Systems|^Genetic Information Processing|^Cellular Processes",
           pathway_class),
    !is.na(p_adjust),
    p_adjust > 1e-3
  )


# Generate pathway error bar plot
picrust_600 <- pathway_errorbar(
  abundance = kegg_abundance_picrust_600,
  daa_results_df = daa_600_filtered,
  Group = metadata_picrust_600$Sample_Type,
  p_values_threshold = 0.05,
  order = "pathway_class",
  select = NULL,
  ko_to_kegg = TRUE,
  p_value_bar = TRUE,
  colors = sample_colors,
  x_lab = "pathway_name"
)

picrust_600 

#save as image
ggsave("picrust_600 .png", plot = picrust_600, 
       path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 45, height = 15, units = c("in"), device = "png")

#save as pdf
ggsave("picrust_600 .pdf", plot = picrust_600, 
       path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 45, height = 15, units = c("in"), device = "pdf")


# Save DAA results
write_csv(
  daa_annotated_sub_method_results_df_picrust_600,
  file.path("C:/Aduragbemi/Manuscript/Review 5/picrust/output/daa_annotated_results_DESeq2_picrust_600.csv")
)



#For GDD_1400

#read metadata
metadata_picrust_1400 <- read.csv("C:/Aduragbemi/Manuscript/Review 5/picrust/input/picrust_1400/all_picrust_1400_metadata.csv")

# Make Growing_Degree_Days a factor
metadata_picrust_1400$Sample_Type <- as.factor(metadata_picrust_1400$Sample_Type)

# Verify the grouping variable
print(table(metadata_picrust_1400$Sample_Type))

# Convert KO abundance to KEGG pathway abundance
kegg_abundance_picrust_1400 <- ko2kegg_abundance("C:/Aduragbemi/Manuscript/Review 5/picrust/input/picrust_1400/all_picrust_1400.tsv") 

# Perform pathway DAA using DESeq2
daa_results_df_picrust_1400 <- pathway_daa(
  abundance = kegg_abundance_picrust_1400, 
  metadata = metadata_picrust_1400, 
  group = "Sample_Type", 
  daa_method = "DESeq2", select = NULL, 
  reference = NULL)

# Filter for DESeq2 results
daa_sub_method_results_df_picrust_1400 <- daa_results_df_picrust_1400[daa_results_df_picrust_1400$method == "DESeq2", ]

# Annotate pathway results
daa_annotated_sub_method_results_df_picrust_1400 <- pathway_annotation(
  pathway = "KO", daa_results_df = 
    daa_sub_method_results_df_picrust_1400, 
  ko_to_kegg = TRUE
)

#Keep only significant KO IDs (adjusted p > 1e-3) ie highly significant
daa_1400_filtered <- daa_annotated_sub_method_results_df_picrust_1400 %>%
  dplyr::filter(
    !grepl("^Human Diseases|^Organismal Systems|^Genetic Information Processing|^Cellular Processes",
           pathway_class),
    !is.na(p_adjust),
    p_adjust > 1e-3
)


# Generate pathway error bar plot
picrust_1400 <- pathway_errorbar(
  abundance = kegg_abundance_picrust_1400,
  daa_results_df = daa_1400_filtered,
  Group = metadata_picrust_1400$Sample_Type,
  p_values_threshold = 0.05,
  order = "pathway_class",
  select = NULL,
  ko_to_kegg = TRUE,
  p_value_bar = TRUE,
  colors = sample_colors,
  x_lab = "pathway_name"
)

picrust_1400  

#save as image
ggsave("picrust_1400 .png", plot = picrust_1400, 
       path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 45, height = 15, units = c("in"), device = "png")

#save as pdf
ggsave("picrust_1400 .pdf", plot = picrust_1400, 
       path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 45, height = 15, units = c("in"), device = "pdf")


# Save DAA results
write_csv(
  daa_annotated_sub_method_results_df_picrust_1400,
  file.path("C:/Aduragbemi/Manuscript/Review 5/picrust/output/daa_annotated_results_DESeq2_picrust_1400.csv")
)

