Georgia Microbiome

This is the workflow for Georgia Microbiome analysis across various field location in Georgia. It contain the data output after processing raw data on QIIME2 and build ASV table.

File input

1.Data file input: this script imports raw microbiome data files (metadata, ASV table, taxonomy, FASTA sequences, and phylogenetic tree) into R and formats them for downstream analysis using phyloseq.

Preprocessing

2.Processing microbiome data: This script constructs a clean and decontaminated phyloseq object from processed microbiome RDS files, applies taxonomic and read-depth filtering, and removes contaminants. It then subsets data by tissue type, evaluates sequencing depth distributions, and exports finalized non-normalized phyloseq objects for downstream diversity, network, and differential abundance analyses.

Map construction

3. Map construction for sampling sites: This script generates a publication-quality geographic map of maize sampling sites using ggmap with Stadia Maps terrain basemap tiles. It imports GPS coordinates, defines spatial extent, overlays color-coded sampling locations, customizes legend labels and aesthetics.

Abundance plot

4.Relative Abundance: This script calculates and visualizes bacterial community composition at the phylum level across sample types and growing degree days using cleaned phyloseq objects.

Diversity metrics 

5a. Alpha diversity for sample types
5b. Alpha diversity Location
These scripts performs rarefaction and calculates alpha diversity metrics (Observed ASVs, Shannon diversity, and Faith’s Phylogenetic Diversity) across sample types, growing degree days, and across field locations. It evaluates normality, applies Kruskal–Wallis and Dunn’s post hoc tests, and generates publication-ready violin/boxplots for comparative diversity analyses.
5c. Alpha diversity for environmental factors, soil chemistry and geography: This script rarefies phyloseq datasets and evaluates alpha diversity (Observed richness, Shannon, and Faith’s PD) across soil, rhizosphere, root, and stalk samples. It applies generalized linear models (GLMs) to test relationships between microbial diversity and environmental variables (climate, geography) as well as detailed soil physicochemical properties.
6. Beta Diversity: This script computes beta diversity using weighted UniFrac distances and performs PCoA ordination across all sample types and within each tissue type. It visualizes microbial community differences by sample type and growing degree days (GDD).
7. Geographic distance-unifrac correlation: This script tests distance–decay relationships by correlating geographic distance (Haversine) with weighted UniFrac beta diversity across soil, rhizosphere, root, and stalk microbiomes. It computes Spearman correlations, visualizes regression trends.


Assessment of Noise in Microbiome Sampling

8a. Noise_assessment_functions
8b. Noise_assessment_soil
8c. Noise_assessment_rhizosphere
8d. Noise_assessment_root
8e. Noise_assessment_stalk
8f. Generate_all_figures for noise assessment
These scripts defines a universal noise assessment framework for microbiome data, evaluating alpha and beta diversity variability within complete location–GDD sample pairs. It performs rarefaction, computes diversity-based noise metrics (including within-pair UniFrac distances), compares against random expectations, generates PCoA coordinates, and exports comprehensive summary outputs for reproducible quality assessment.

Environmental descriptors and their association with the microbiome communities

9. Canonical Correspondence Analysis: This script performs CCA-based ordination to evaluate the influence of environmental variables, soil chemistry, and geographic factors on soil, rhizosphere, and root microbial communities.	
10. Mantel Test: This script assesses the relationship between microbial community dissimilarity (weighted UniFrac) and geography, environmental, and soil chemical distance matrices using permutation-based Spearman Mantel tests across soil, rhizosphere, root, and stalk compartments.	
11. Pearson Correlations between the relative abundance of Individual taxa with environmental factors: This script quantifies pearson correlations between bacterial class-level relative abundances and environmental and soil chemistry variables across soil, rhizosphere, root, and stalk compartments.

Differential abundance analysis of taxa 

12a. Differential abundance analysis-class level
12b. Differential abundance analysis - volcano plot genus level
These scripts performs differential abundance analysis using DESeq2 to identify bacterial classes and genus significantly enriched between 600 and 1400 Growing Degree Days (GDD) across soil, rhizosphere, root, and stalk compartments. The plot were at class and genus level respectively.

Network co-occurrence analysis 

13. Network analysis of sample types: This script constructs and compares co-occurrence networks of soil, rhizosphere, root, and stalk microbiomes at 600 and 1400 Growing Degree Days (GDD) using Spearman correlations.
Metagenomic Prediction of microbial communities
14. Picrust2 metagenome prediction: This script prepares ASV sequences, OTU tables, and metadata for PICRUSt2 functional pathway prediction, then performs differential pathway abundance analysis (DESeq2-based) across Growing Degree Days (600 vs 1400) and among sample types (Soil, Rhizosphere, Root, Stalk)
