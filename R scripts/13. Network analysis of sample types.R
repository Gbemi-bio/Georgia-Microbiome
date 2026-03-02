#Load required libraries
library(file2meco)
library(microeco)
library(phyloseq)
library(rgexf)
library(tidyverse)
library(ggpubr)
library(meconetcomp)

# ================================
# read in the non-normalized data
# ================================
soil_bac_sub_no_bad_Filtered_final_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/soil_bac_sub_no_bad_Filtered_final_non_normalized.rds")
rhizo_bac_sub_no_bad_Filtered_final_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/rhizo_bac_sub_no_bad_Filtered_final_non_normalized.rds")
root_bac_sub_no_bad_Filtered_final_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/root_bac_sub_no_bad_Filtered_final_non_normalized.rds")
stalk_bac_sub_no_bad_Filtered_final_non_normalized <- readRDS(file = "C:/Aduragbemi/Microbiome/RDS/stalk_bac_sub_no_bad_Filtered_final_non_normalized.rds")


# ================================
# Network analysis of sample types
# ================================

#For Soil samples

#Normalized counts to relative abundance and filtered low-abundance taxa
soil_norm = transform_sample_counts(soil_bac_sub_no_bad_Filtered_final_non_normalized, function(x) x / sum(x) )
soil_network = filter_taxa(soil_norm, function(x) mean(x) > 0.0001, TRUE)
soil_network


#Growing degree 600
#subset the soil sample to 600 growing degree days
soil_network_600 <- subset_samples(soil_network, Growing_Degree_Days == "600")
soil_network_600


#first convert to microtable
soil_amp_600 <- phyloseq2meco(soil_network_600)
soil_amp_600

#Construct Spearman correlation network using WGCNA
net_soil_600 <- trans_network$new(dataset = soil_amp_600, 
                                  cor_method = "spearman", 
                                  use_WGCNA_pearson_spearman = TRUE, 
                                  filter_thres = 0.0001)
net_soil_600$cal_network(COR_p_thres = 0.05, COR_cut = 0.6, COR_p_adjust = "fdr")


#Identify network structure using fast greedy community detection algorithm
net_soil_600$cal_module(method = "cluster_fast_greedy")

#Saved constructed soil microbial network as GEXF file
net_soil_600$save_network(filepath = "C:/Aduragbemi/Manuscript/Review 5/Network files/soil_network_600.gexf")

# calculate network attributes and export results to CSV.
net_soil_attr_600 <- net_soil_600$cal_network_attr()
net_soil_attr_600 <- net_soil_600$res_network_attr
write.csv(net_soil_attr_600, "C:/Aduragbemi/Manuscript/Review 5/Network files/net_soil_attr_600.csv")


# get node properties
net_soil_600$get_node_table(node_roles = TRUE)
net_soil_600$res_node_table

#Extract node-level results table from constructed network object and write to csv
soil_res_node_table_600 <- net_soil_600[["res_node_table"]]
write.csv(soil_res_node_table_600, "C:/Aduragbemi/Manuscript/Review 5/Network files/soil_res_node_table_600.csv")

# get edge properties
net_soil_600$get_edge_table()
net_soil_600$res_edge_table 

#get network adjacency matrix
net_soil_600$get_adjacency_matrix()
net_soil_600$res_adjacency_matrix


#Visualize taxa roles with customized connectivity thresholds and publication styling
soil_connectivity_600 <- net_soil_600$plot_taxa_roles(
  use_type = 1,
  roles_color_background = FALSE,
  size = 9
)

soil_connectivity_600 <- soil_connectivity_600 +
  theme(
    plot.title = element_text(
      family = "serif",
      face = "bold",
      size = 35
    ),
    axis.title = element_text(
      family = "serif",
      size = 45,
      colour = "black",
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      size = 45,
      colour = "black"
    ),
    axis.text.y = element_text(
      angle = 0,
      hjust = 1,
      vjust = 0.5,
      size = 45,
      colour = "black"
    ),
    legend.text = element_text(size = 45),
    legend.title = element_text(
      face = "bold",
      size = 45
    ),
    legend.key.size = unit(2, "cm"),
    strip.text = element_text(
      size = 45,
      face = "bold",
      angle = 0
    ),
    axis.line = element_line(
      colour = "black",
      size = 1,
      linetype = "solid"
    )
  ) +
  geom_vline(
    xintercept = 0.62,
    linetype = "solid",
    color = "black",
    size = 1
  ) +
  geom_hline(
    yintercept = 2.5,
    linetype = "solid",
    color = "black",
    size = 1
  )

soil_connectivity_600


#save plot
ggsave("soil_connectivity_600.png", plot = soil_connectivity_600 , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "png")

ggsave("soil_connectivity_600.pdf", plot = soil_connectivity_600 , path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "pdf")





#For Growing Degree 1400
#subset the soil sample to 1400 growing degree days
soil_network_1400 <- subset_samples(soil_network, Growing_Degree_Days == "1400")
soil_network_1400


#first convert to microtable
soil_amp_1400 <- phyloseq2meco(soil_network_1400)
soil_amp_1400

#Construct Spearman correlation network using WGCNA
net_soil_1400 <- trans_network$new(dataset = soil_amp_1400, 
                                   cor_method = "spearman", 
                                   use_WGCNA_pearson_spearman = TRUE, 
                                   filter_thres = 0.0001)
net_soil_1400$cal_network(COR_p_thres = 0.05, COR_cut = 0.6, COR_p_adjust = "fdr")


#Identify network structure using fast greedy community detection algorithm
net_soil_1400$cal_module(method = "cluster_fast_greedy")

#Save constructed soil microbial network as GEXF file
net_soil_1400$save_network(filepath = "C:/Aduragbemi/Manuscript/Review 5/Network files/soil_network_1400.gexf")

# calculate network attributes and export results to CSV.
net_soil_attr_1400 <- net_soil_1400$cal_network_attr()
net_soil_attr_1400 <- net_soil_1400$res_network_attr
write.csv(net_soil_attr_1400, "C:/Aduragbemi/Manuscript/Review 5/Network files/net_soil_attr_1400.csv")

# get node properties
net_soil_1400$get_node_table(node_roles = TRUE)
net_soil_1400$res_node_table

#Extract node-level results table from constructed network object and write to csv
soil_res_node_table_1400 <- net_soil_1400[["res_node_table"]]
write.csv(soil_res_node_table_1400, "C:/Aduragbemi/Manuscript/Review 5/Network files/soil_res_node_table_1400.csv")

# get edge properties
net_soil_1400$get_edge_table()
net_soil_1400$res_edge_table 

#get network adjacency matrix
net_soil_1400$get_adjacency_matrix()
net_soil_1400$res_adjacency_matrix


#Visualized taxa roles with customized connectivity thresholds and publication styling
soil_connectivity_1400 <- net_soil_1400$plot_taxa_roles(
  use_type = 1,
  roles_color_background = FALSE,
  size = 9
)

soil_connectivity_1400 <- soil_connectivity_1400 +
  theme(
    plot.title = element_text(
      family = "serif",
      face = "bold",
      size = 35
    ),
    axis.title = element_text(
      family = "serif",
      size = 45,
      colour = "black",
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      size = 45,
      colour = "black"
    ),
    axis.text.y = element_text(
      angle = 0,
      hjust = 1,
      vjust = 0.5,
      size = 45,
      colour = "black"
    ),
    legend.text = element_text(size = 45),
    legend.title = element_text(
      face = "bold",
      size = 45
    ),
    legend.key.size = unit(2, "cm"),
    strip.text = element_text(
      size = 45,
      face = "bold",
      angle = 0
    ),
    axis.line = element_line(
      colour = "black",
      size = 1,
      linetype = "solid"
    )
  ) +
  geom_vline(
    xintercept = 0.62,
    linetype = "solid",
    color = "black",
    size = 1
  ) +
  geom_hline(
    yintercept = 2.5,
    linetype = "solid",
    color = "black",
    size = 1
  )

soil_connectivity_1400


#save plot
ggsave("soil_connectivity_1400.png", plot = soil_connectivity_1400 , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "png")

ggsave("soil_connectivity_1400.pdf", plot = soil_connectivity_1400 , path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "pdf")



#For Rhizosphere samples

#Normalized counts to relative abundance and filtered low-abundance taxa
rhizo_norm = transform_sample_counts(rhizo_bac_sub_no_bad_Filtered_final_non_normalized, function(x) x / sum(x) )
rhizo_network = filter_taxa(rhizo_norm, function(x) mean(x) > 0.0001, TRUE)
rhizo_network


#Growing degree 600
#subset the rhizo sample to 600 growing degree days
rhizo_network_600 <- subset_samples(rhizo_network, Growing_Degree_Days == "600")
rhizo_network_600


#first convert to microtable
rhizo_amp_600 <- phyloseq2meco(rhizo_network_600)
rhizo_amp_600

#Construct Spearman correlation network using WGCNA
net_rhizo_600 <- trans_network$new(dataset = rhizo_amp_600, 
                                  cor_method = "spearman", 
                                  use_WGCNA_pearson_spearman = TRUE, 
                                  filter_thres = 0.0001)
net_rhizo_600$cal_network(COR_p_thres = 0.05, COR_cut = 0.6, COR_p_adjust = "fdr")


#Identify network structure using fast greedy community detection algorithm
net_rhizo_600$cal_module(method = "cluster_fast_greedy")

#Saved constructed rhizo microbial network as GEXF file
net_rhizo_600$save_network(filepath = "C:/Aduragbemi/Manuscript/Review 5/Network files/rhizo_network_600.gexf")

# calculate network attributes and export results to CSV.
net_rhizo_attr_600 <- net_rhizo_600$cal_network_attr()
net_rhizo_attr_600 <- net_rhizo_600$res_network_attr
write.csv(net_rhizo_attr_600, "C:/Aduragbemi/Manuscript/Review 5/Network files/net_rhizo_attr_600.csv")


# get node properties
net_rhizo_600$get_node_table(node_roles = TRUE)
net_rhizo_600$res_node_table

#Extract node-level results table from constructed network object and write to csv
rhizo_res_node_table_600 <- net_rhizo_600[["res_node_table"]]
write.csv(rhizo_res_node_table_600, "C:/Aduragbemi/Manuscript/Review 5/Network files/rhizo_res_node_table_600.csv")

# get edge properties
net_rhizo_600$get_edge_table()
net_rhizo_600$res_edge_table 

#get network adjacency matrix
net_rhizo_600$get_adjacency_matrix()
net_rhizo_600$res_adjacency_matrix


#Visualize taxa roles with customized connectivity thresholds and publication styling
rhizo_connectivity_600 <- net_rhizo_600$plot_taxa_roles(
  use_type = 1,
  roles_color_background = FALSE,
  size = 9
)

rhizo_connectivity_600 <- rhizo_connectivity_600 +
  theme(
    plot.title = element_text(
      family = "serif",
      face = "bold",
      size = 35
    ),
    axis.title = element_text(
      family = "serif",
      size = 45,
      colour = "black",
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      size = 45,
      colour = "black"
    ),
    axis.text.y = element_text(
      angle = 0,
      hjust = 1,
      vjust = 0.5,
      size = 45,
      colour = "black"
    ),
    legend.text = element_text(size = 45),
    legend.title = element_text(
      face = "bold",
      size = 45
    ),
    legend.key.size = unit(2, "cm"),
    strip.text = element_text(
      size = 45,
      face = "bold",
      angle = 0
    ),
    axis.line = element_line(
      colour = "black",
      size = 1,
      linetype = "solid"
    )
  ) +
  geom_vline(
    xintercept = 0.62,
    linetype = "solid",
    color = "black",
    size = 1
  ) +
  geom_hline(
    yintercept = 2.5,
    linetype = "solid",
    color = "black",
    size = 1
  )

rhizo_connectivity_600


#save plot
ggsave("rhizo_connectivity_600.png", plot = rhizo_connectivity_600 , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "png")

ggsave("rhizo_connectivity_600.pdf", plot = rhizo_connectivity_600 , path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "pdf")





#For Growing Degree 1400
#subset the rhizo sample to 1400 growing degree days
rhizo_network_1400 <- subset_samples(rhizo_network, Growing_Degree_Days == "1400")
rhizo_network_1400


#first convert to microtable
rhizo_amp_1400 <- phyloseq2meco(rhizo_network_1400)
rhizo_amp_1400

#Construct Spearman correlation network using WGCNA
net_rhizo_1400 <- trans_network$new(dataset = rhizo_amp_1400, 
                                   cor_method = "spearman", 
                                   use_WGCNA_pearson_spearman = TRUE, 
                                   filter_thres = 0.0001)
net_rhizo_1400$cal_network(COR_p_thres = 0.05, COR_cut = 0.6, COR_p_adjust = "fdr")


#Identify network structure using fast greedy community detection algorithm
net_rhizo_1400$cal_module(method = "cluster_fast_greedy")

#Save constructed rhizo microbial network as GEXF file
net_rhizo_1400$save_network(filepath = "C:/Aduragbemi/Manuscript/Review 5/Network files/rhizo_network_1400.gexf")

# calculate network attributes and export results to CSV.
net_rhizo_attr_1400 <- net_rhizo_1400$cal_network_attr()
net_rhizo_attr_1400 <- net_rhizo_1400$res_network_attr
write.csv(net_rhizo_attr_1400, "C:/Aduragbemi/Manuscript/Review 5/Network files/net_rhizo_attr_1400.csv")

# get node properties
net_rhizo_1400$get_node_table(node_roles = TRUE)
net_rhizo_1400$res_node_table

#Extract node-level results table from constructed network object and write to csv
rhizo_res_node_table_1400 <- net_rhizo_1400[["res_node_table"]]
write.csv(rhizo_res_node_table_1400, "C:/Aduragbemi/Manuscript/Review 5/Network files/rhizo_res_node_table_1400.csv")

# get edge properties
net_rhizo_1400$get_edge_table()
net_rhizo_1400$res_edge_table 

#get network adjacency matrix
net_rhizo_1400$get_adjacency_matrix()
net_rhizo_1400$res_adjacency_matrix


#Visualized taxa roles with customized connectivity thresholds and publication styling
rhizo_connectivity_1400 <- net_rhizo_1400$plot_taxa_roles(
  use_type = 1,
  roles_color_background = FALSE,
  size = 9
)

rhizo_connectivity_1400 <- rhizo_connectivity_1400 +
  theme(
    plot.title = element_text(
      family = "serif",
      face = "bold",
      size = 35
    ),
    axis.title = element_text(
      family = "serif",
      size = 45,
      colour = "black",
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      size = 45,
      colour = "black"
    ),
    axis.text.y = element_text(
      angle = 0,
      hjust = 1,
      vjust = 0.5,
      size = 45,
      colour = "black"
    ),
    legend.text = element_text(size = 45),
    legend.title = element_text(
      face = "bold",
      size = 45
    ),
    legend.key.size = unit(2, "cm"),
    strip.text = element_text(
      size = 45,
      face = "bold",
      angle = 0
    ),
    axis.line = element_line(
      colour = "black",
      size = 1,
      linetype = "solid"
    )
  ) +
  geom_vline(
    xintercept = 0.62,
    linetype = "solid",
    color = "black",
    size = 1
  ) +
  geom_hline(
    yintercept = 2.5,
    linetype = "solid",
    color = "black",
    size = 1
  )

rhizo_connectivity_1400


#save plot
ggsave("rhizo_connectivity_1400.png", plot = rhizo_connectivity_1400 , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "png")

ggsave("rhizo_connectivity_1400.pdf", plot = rhizo_connectivity_1400 , path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "pdf")




#For Root samples

#Normalized counts to relative abundance and filtered low-abundance taxa
root_norm = transform_sample_counts(root_bac_sub_no_bad_Filtered_final_non_normalized, function(x) x / sum(x) )
root_network = filter_taxa(root_norm, function(x) mean(x) > 0.0001, TRUE)
root_network


#Growing degree 600
#subset the root sample to 600 growing degree days
root_network_600 <- subset_samples(root_network, Growing_Degree_Days == "600")
root_network_600


#first convert to microtable
root_amp_600 <- phyloseq2meco(root_network_600)
root_amp_600

#Construct Spearman correlation network using WGCNA
net_root_600 <- trans_network$new(dataset = root_amp_600, 
                                   cor_method = "spearman", 
                                   use_WGCNA_pearson_spearman = TRUE, 
                                   filter_thres = 0.0001)
net_root_600$cal_network(COR_p_thres = 0.05, COR_cut = 0.6, COR_p_adjust = "fdr")


#Identify network structure using fast greedy community detection algorithm
net_root_600$cal_module(method = "cluster_fast_greedy")

#Saved constructed root microbial network as GEXF file
net_root_600$save_network(filepath = "C:/Aduragbemi/Manuscript/Review 5/Network files/root_network_600.gexf")

# calculate network attributes and export results to CSV.
net_root_attr_600 <- net_root_600$cal_network_attr()
net_root_attr_600 <- net_root_600$res_network_attr
write.csv(net_root_attr_600, "C:/Aduragbemi/Manuscript/Review 5/Network files/net_root_attr_600.csv")


# get node properties
net_root_600$get_node_table(node_roles = TRUE)
net_root_600$res_node_table

#Extract node-level results table from constructed network object and write to csv
root_res_node_table_600 <- net_root_600[["res_node_table"]]
write.csv(root_res_node_table_600, "C:/Aduragbemi/Manuscript/Review 5/Network files/root_res_node_table_600.csv")

# get edge properties
net_root_600$get_edge_table()
net_root_600$res_edge_table 

#get network adjacency matrix
net_root_600$get_adjacency_matrix()
net_root_600$res_adjacency_matrix


#Visualize taxa roles with customized connectivity thresholds and publication styling
root_connectivity_600 <- net_root_600$plot_taxa_roles(
  use_type = 1,
  roles_color_background = FALSE,
  size = 9
)

root_connectivity_600 <- root_connectivity_600 +
  theme(
    plot.title = element_text(
      family = "serif",
      face = "bold",
      size = 35
    ),
    axis.title = element_text(
      family = "serif",
      size = 45,
      colour = "black",
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      size = 45,
      colour = "black"
    ),
    axis.text.y = element_text(
      angle = 0,
      hjust = 1,
      vjust = 0.5,
      size = 45,
      colour = "black"
    ),
    legend.text = element_text(size = 45),
    legend.title = element_text(
      face = "bold",
      size = 45
    ),
    legend.key.size = unit(2, "cm"),
    strip.text = element_text(
      size = 45,
      face = "bold",
      angle = 0
    ),
    axis.line = element_line(
      colour = "black",
      size = 1,
      linetype = "solid"
    )
  ) +
  geom_vline(
    xintercept = 0.62,
    linetype = "solid",
    color = "black",
    size = 1
  ) +
  geom_hline(
    yintercept = 2.5,
    linetype = "solid",
    color = "black",
    size = 1
  )

root_connectivity_600


#save plot
ggsave("root_connectivity_600.png", plot = root_connectivity_600 , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "png")

ggsave("root_connectivity_600.pdf", plot = root_connectivity_600 , path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "pdf")





#For Growing Degree 1400
#subset the root sample to 1400 growing degree days
root_network_1400 <- subset_samples(root_network, Growing_Degree_Days == "1400")
root_network_1400


#first convert to microtable
root_amp_1400 <- phyloseq2meco(root_network_1400)
root_amp_1400

#Construct Spearman correlation network using WGCNA
net_root_1400 <- trans_network$new(dataset = root_amp_1400, 
                                    cor_method = "spearman", 
                                    use_WGCNA_pearson_spearman = TRUE, 
                                    filter_thres = 0.0001)
net_root_1400$cal_network(COR_p_thres = 0.05, COR_cut = 0.6, COR_p_adjust = "fdr")


#Identify network structure using fast greedy community detection algorithm
net_root_1400$cal_module(method = "cluster_fast_greedy")

#Save constructed root microbial network as GEXF file
net_root_1400$save_network(filepath = "C:/Aduragbemi/Manuscript/Review 5/Network files/root_network_1400.gexf")

# calculate network attributes and export results to CSV.
net_root_attr_1400 <- net_root_1400$cal_network_attr()
net_root_attr_1400 <- net_root_1400$res_network_attr
write.csv(net_root_attr_1400, "C:/Aduragbemi/Manuscript/Review 5/Network files/net_root_attr_1400.csv")

# get node properties
net_root_1400$get_node_table(node_roles = TRUE)
net_root_1400$res_node_table

#Extract node-level results table from constructed network object and write to csv
root_res_node_table_1400 <- net_root_1400[["res_node_table"]]
write.csv(root_res_node_table_1400, "C:/Aduragbemi/Manuscript/Review 5/Network files/root_res_node_table_1400.csv")

# get edge properties
net_root_1400$get_edge_table()
net_root_1400$res_edge_table 

#get network adjacency matrix
net_root_1400$get_adjacency_matrix()
net_root_1400$res_adjacency_matrix


#Visualized taxa roles with customized connectivity thresholds and publication styling
root_connectivity_1400 <- net_root_1400$plot_taxa_roles(
  use_type = 1,
  roles_color_background = FALSE,
  size = 9
)

root_connectivity_1400 <- root_connectivity_1400 +
  theme(
    plot.title = element_text(
      family = "serif",
      face = "bold",
      size = 35
    ),
    axis.title = element_text(
      family = "serif",
      size = 45,
      colour = "black",
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      size = 45,
      colour = "black"
    ),
    axis.text.y = element_text(
      angle = 0,
      hjust = 1,
      vjust = 0.5,
      size = 45,
      colour = "black"
    ),
    legend.text = element_text(size = 45),
    legend.title = element_text(
      face = "bold",
      size = 45
    ),
    legend.key.size = unit(2, "cm"),
    strip.text = element_text(
      size = 45,
      face = "bold",
      angle = 0
    ),
    axis.line = element_line(
      colour = "black",
      size = 1,
      linetype = "solid"
    )
  ) +
  geom_vline(
    xintercept = 0.62,
    linetype = "solid",
    color = "black",
    size = 1
  ) +
  geom_hline(
    yintercept = 2.5,
    linetype = "solid",
    color = "black",
    size = 1
  )

root_connectivity_1400


#save plot
ggsave("root_connectivity_1400.png", plot = root_connectivity_1400 , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "png")

ggsave("root_connectivity_1400.pdf", plot = root_connectivity_1400 , path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "pdf")



#For Stalk samples

#Normalized counts to relative abundance and filtered low-abundance taxa
stalk_norm = transform_sample_counts(stalk_bac_sub_no_bad_Filtered_final_non_normalized, function(x) x / sum(x) )
stalk_network = filter_taxa(stalk_norm, function(x) mean(x) > 0.0001, TRUE)
stalk_network


#Growing degree 600
#subset the stalk sample to 600 growing degree days
stalk_network_600 <- subset_samples(stalk_network, Growing_Degree_Days == "600")
stalk_network_600


#first convert to microtable
stalk_amp_600 <- phyloseq2meco(stalk_network_600)
stalk_amp_600

#Construct Spearman correlation network using WGCNA
net_stalk_600 <- trans_network$new(dataset = stalk_amp_600, 
                                  cor_method = "spearman", 
                                  use_WGCNA_pearson_spearman = TRUE, 
                                  filter_thres = 0.0001)
net_stalk_600$cal_network(COR_p_thres = 0.05, COR_cut = 0.6, COR_p_adjust = "fdr")


#Identify network structure using fast greedy community detection algorithm
net_stalk_600$cal_module(method = "cluster_fast_greedy")

#Saved constructed stalk microbial network as GEXF file
net_stalk_600$save_network(filepath = "C:/Aduragbemi/Manuscript/Review 5/Network files/stalk_network_600.gexf")

# calculate network attributes and export results to CSV.
net_stalk_attr_600 <- net_stalk_600$cal_network_attr()
net_stalk_attr_600 <- net_stalk_600$res_network_attr
write.csv(net_stalk_attr_600, "C:/Aduragbemi/Manuscript/Review 5/Network files/net_stalk_attr_600.csv")


# get node properties
net_stalk_600$get_node_table(node_roles = TRUE)
net_stalk_600$res_node_table

#Extract node-level results table from constructed network object and write to csv
stalk_res_node_table_600 <- net_stalk_600[["res_node_table"]]
write.csv(stalk_res_node_table_600, "C:/Aduragbemi/Manuscript/Review 5/Network files/stalk_res_node_table_600.csv")

# get edge properties
net_stalk_600$get_edge_table()
net_stalk_600$res_edge_table 

#get network adjacency matrix
net_stalk_600$get_adjacency_matrix()
net_stalk_600$res_adjacency_matrix


#Visualize taxa roles with customized connectivity thresholds and publication styling
stalk_connectivity_600 <- net_stalk_600$plot_taxa_roles(
  use_type = 1,
  roles_color_background = FALSE,
  size = 9
)

stalk_connectivity_600 <- stalk_connectivity_600 +
  theme(
    plot.title = element_text(
      family = "serif",
      face = "bold",
      size = 35
    ),
    axis.title = element_text(
      family = "serif",
      size = 45,
      colour = "black",
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      size = 45,
      colour = "black"
    ),
    axis.text.y = element_text(
      angle = 0,
      hjust = 1,
      vjust = 0.5,
      size = 45,
      colour = "black"
    ),
    legend.text = element_text(size = 45),
    legend.title = element_text(
      face = "bold",
      size = 45
    ),
    legend.key.size = unit(2, "cm"),
    strip.text = element_text(
      size = 45,
      face = "bold",
      angle = 0
    ),
    axis.line = element_line(
      colour = "black",
      size = 1,
      linetype = "solid"
    )
  ) +
  geom_vline(
    xintercept = 0.62,
    linetype = "solid",
    color = "black",
    size = 1
  ) +
  geom_hline(
    yintercept = 2.5,
    linetype = "solid",
    color = "black",
    size = 1
  )

stalk_connectivity_600


#save plot
ggsave("stalk_connectivity_600.png", plot = stalk_connectivity_600 , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "png")

ggsave("stalk_connectivity_600.pdf", plot = stalk_connectivity_600 , path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "pdf")





#For Growing Degree 1400
#subset the stalk sample to 1400 growing degree days
stalk_network_1400 <- subset_samples(stalk_network, Growing_Degree_Days == "1400")
stalk_network_1400


#first convert to microtable
stalk_amp_1400 <- phyloseq2meco(stalk_network_1400)
stalk_amp_1400

#Construct Spearman correlation network using WGCNA
net_stalk_1400 <- trans_network$new(dataset = stalk_amp_1400, 
                                   cor_method = "spearman", 
                                   use_WGCNA_pearson_spearman = TRUE, 
                                   filter_thres = 0.0001)
net_stalk_1400$cal_network(COR_p_thres = 0.05, COR_cut = 0.6, COR_p_adjust = "fdr")


#Identify network structure using fast greedy community detection algorithm
net_stalk_1400$cal_module(method = "cluster_fast_greedy")

#Save constructed stalk microbial network as GEXF file
net_stalk_1400$save_network(filepath = "C:/Aduragbemi/Manuscript/Review 5/Network files/stalk_network_1400.gexf")

# calculate network attributes and export results to CSV.
net_stalk_attr_1400 <- net_stalk_1400$cal_network_attr()
net_stalk_attr_1400 <- net_stalk_1400$res_network_attr
write.csv(net_stalk_attr_1400, "C:/Aduragbemi/Manuscript/Review 5/Network files/net_stalk_attr_1400.csv")

# get node properties
net_stalk_1400$get_node_table(node_roles = TRUE)
net_stalk_1400$res_node_table

#Extract node-level results table from constructed network object and write to csv
stalk_res_node_table_1400 <- net_stalk_1400[["res_node_table"]]
write.csv(stalk_res_node_table_1400, "C:/Aduragbemi/Manuscript/Review 5/Network files/stalk_res_node_table_1400.csv")

# get edge properties
net_stalk_1400$get_edge_table()
net_stalk_1400$res_edge_table 

#get network adjacency matrix
net_stalk_1400$get_adjacency_matrix()
net_stalk_1400$res_adjacency_matrix


#Visualized taxa roles with customized connectivity thresholds and publication styling
stalk_connectivity_1400 <- net_stalk_1400$plot_taxa_roles(
  use_type = 1,
  roles_color_background = FALSE,
  size = 9
)

stalk_connectivity_1400 <- stalk_connectivity_1400 +
  theme(
    plot.title = element_text(
      family = "serif",
      face = "bold",
      size = 35
    ),
    axis.title = element_text(
      family = "serif",
      size = 45,
      colour = "black",
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      size = 45,
      colour = "black"
    ),
    axis.text.y = element_text(
      angle = 0,
      hjust = 1,
      vjust = 0.5,
      size = 45,
      colour = "black"
    ),
    legend.text = element_text(size = 45),
    legend.title = element_text(
      face = "bold",
      size = 45
    ),
    legend.key.size = unit(2, "cm"),
    strip.text = element_text(
      size = 45,
      face = "bold",
      angle = 0
    ),
    axis.line = element_line(
      colour = "black",
      size = 1,
      linetype = "solid"
    )
  ) +
  geom_vline(
    xintercept = 0.62,
    linetype = "solid",
    color = "black",
    size = 1
  ) +
  geom_hline(
    yintercept = 2.5,
    linetype = "solid",
    color = "black",
    size = 1
  )

stalk_connectivity_1400


#save plot
ggsave("stalk_connectivity_1400.png", plot = stalk_connectivity_1400 , path = "C:/Aduragbemi/Manuscript/Review 5/IMAGES5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "png")

ggsave("stalk_connectivity_1400.pdf", plot = stalk_connectivity_1400 , path = "C:/Aduragbemi/Manuscript/Review 5/PDF5", dpi = 700, 
       width = 15, height = 10, units = c("in"), device = "pdf")

