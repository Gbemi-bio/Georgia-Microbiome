#Load Libraries#
library(phyloseq)

# ===========================================================================
# I need to first input our data file to R
# Make a few edits first, then save as RDS before continuing analysis in R
# ===========================================================================

#First the metadata#
samp_dat_bac = read.delim("C:/Aduragbemi/Microbiome/raw files/metadata-use.tsv")
samp_dat_bac
rownames(samp_dat_bac) <- samp_dat_bac$Samples #row names must match OTU table headers
SAMP.bac <- phyloseq::sample_data(samp_dat_bac)
SAMP.bac

#save metadata as RDS#
saveRDS(SAMP.bac, file = "C:/Aduragbemi/Microbiome/RDS/SAMP.bac.rds")
SAMP.bac <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/SAMP.bac.rds")
SAMP.bac

#Next, the OTU Table#
otu_bac <- read.csv("C:/Aduragbemi/Microbiome/raw files/asv-table.csv")
otu_bac
rownames(otu_bac) <- otu_bac$OTU
otu_bac <- otu_bac[,-1]
otu_bac
OTU.bac <- phyloseq::otu_table(otu_bac, taxa_are_rows = TRUE)
OTU.bac
any(is.na(otu_bac)) # no NA in the OTU table

#save OTU Table as RDS#
saveRDS(OTU.bac, file = "C:/Aduragbemi/Microbiome/RDS/OTU.bac.rds")
OTU.bac <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/OTU.bac.rds")
OTU.bac


#Next, the Taxonomy#
taxonomy.bac <- read.csv("C:/Aduragbemi/Microbiome/taxonomy.csv")
head(taxonomy.bac)

#replace the uncultured in the Genus column, with uncultured and phylum name
for (i in 1:nrow(taxonomy.bac)) {
  if (taxonomy.bac$Genus[i] == "uncultured") {
    taxonomy.bac$Genus[i] <- paste0("uncultured ", taxonomy.bac$Phylum[i], sep = " ")
  }
}
head(taxonomy.bac)
rownames(taxonomy.bac) <- taxonomy.bac$OTU
TAX.bac <- phyloseq::tax_table(as.matrix(taxonomy.bac))
TAX.bac
head(TAX.bac)

all.equal(rownames(samp_dat_bac), colnames(otu_bac))

#save taxa as RDS#
saveRDS(TAX.bac, file = "C:/Aduragbemi/Microbiome/RDS/taxonomy.rds")
TAX.bac <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/taxonomy.rds")
TAX.bac

#Next, is Fasta#
FASTA.bac <- readDNAStringSet("C:/Aduragbemi/Microbiome/raw files/dna-sequences.fasta", 
                              format="fasta", seek.first.rec=TRUE, use.names=TRUE)
FASTA.bac

#save fasta as RDS#
saveRDS(FASTA.bac, file = "C:/Aduragbemi/Microbiome/RDS/FASTA.bac.rds")
FASTA.bac <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/FASTA.bac.rds")
FASTA.bac


#Next, is Phylogenetic tree#
tree <- phyloseq::read_tree("C:/Aduragbemi/Microbiome/raw files/rooted_tree.nwk")

#save tree as RDS#
saveRDS(tree, file = "C:/Aduragbemi/Microbiome/RDS/tree.rds")
tree <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/tree.rds")
tree










