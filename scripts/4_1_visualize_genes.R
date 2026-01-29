library(Seurat)
library(tidyverse)
library(sf)
library(spdep)

dataset <- "_caw009" #or "_caw009" or ""
data_path <- paste0("~/visium-hd/data/Visium_HD_Liver", toupper(dataset), "/")
plot_path <- paste0("plots/")


color_palette <- c("Hepatocytes" = "#B4B5B5FF",
                   "CentralVeinEndothelialcells" = "#FED8B1FF",
                   "LSECs" = "#FBB05FFF",
                   "PortalVeinEndothelialcells" = "#CC7722FF",
                   "LymphaticEndothelialcells" = "#8F4716FF",
                   "Cholangiocytes" = "#C61B84FF",
                   "HsPCs" = "#F19FC3FF",
                   "Stellatecells" = "#A31A2AFF",
                   "Mesothelialcells" = "#D0110BFF",
                   "Fibroblasts" = "#E45466FF",
                   "CapsularFibroblasts" = "#D46F6CFF",
                   "Kupffercells" = "#5DA6DBFF",
                   "MonocytesMonocytederivedcells" = "#a3daf3",
                   "cDC1s" = "#893A86FF",
                   "cDC2s" = "#893A86FF",
                   "pDCs" = "#893A86FF",
                   "MigcDCs" = "#893A86FF",
                   "Bcells" = "#9C7EBAFF",
                   "NKcells" = "#4A6E34FF",
                   "Tcells" = "#3AB04AFF",
                   "ILC1s" = "#A3D7BAFF",
                   "Basophils" = "#191919",
                   "Neutrophils" = "#727272")

celltype_order <- names(color_palette)

bin_size <- 8 # 8, 16, or 32

# Fill the text with 0
bin_size_str <- sprintf("%03dum", bin_size)

visium_obj <- readRDS(paste0(data_path, "Visium_HD_Liver", toupper(dataset), "_",
                             bin_size_str, ".rds"))

# Plot histogram of "Gdf2"
genes_oi <- c("Gdf2", "Bmp10", "Acvrl1")
genes_oi %in% rownames(visium_obj)

genes_oi_mat <- GetAssayData(visium_obj, layer="counts")[genes_oi, ]

genes_oi_df <- summary(genes_oi_mat) %>%
  `colnames<-`(c("gene", "cell", "count")) %>% 
  # Sum up counts in cell
  group_by(gene, count) %>% 
  summarise(n = n()) %>% 
  mutate(gene = genes_oi[gene]) %>% 
  ungroup() %>% 
  # Add count=0
  complete(gene, count = 0) %>% 
  group_by(gene) %>%
  mutate(n = ifelse(is.na(n), ncol(visium_obj) - sum(n, na.rm = TRUE), n))

genes_oi_df %>% 
  ggplot(aes(x = count, y=n)) +
  geom_bar(stat="identity") +
  facet_wrap(~gene, scales = "free")

# How many bins have both Gdf2 and Acvrl1 > 0?
gdf_acvrl1_df <-  GetAssayData(visium_obj, layer="counts")[c("Gdf2", "Acvrl1"), ] %>% 
  summary() %>% `colnames<-`(c("gene", "cell", "count")) %>% 
  group_by(cell) %>% 
  summarise(genes_expressed = n())

gdf_acvrl1_df$genes_expressed %>% table

# How many bins have both Bmp10 and Acvrl1 > 0?
bmp_acvrl1_df <-  GetAssayData(visium_obj, layer="counts")[c("Bmp10", "Acvrl1"), ] %>% 
  summary() %>% `colnames<-`(c("gene", "cell", "count")) %>% 
  group_by(cell) %>% 
  summarise(genes_expressed = n())

bmp_acvrl1_df$genes_expressed %>% table


# Create null distribution with 1000 random gene pairs
square_size <- visium_obj@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot
set.seed(123)
all_genes <- rownames(visium_obj)
n_random <- 1000
random_gene_pairs <- tibble(
  gene1 = sample(all_genes, n_random, replace = TRUE),
  gene2 = sample(all_genes, n_random, replace = TRUE)
)
all(!random_gene_pairs$gene1 == random_gene_pairs$gene2)

# Get all info 
tissue_coords <- GetTissueCoordinates(visium_obj)
random_genes_mat <- GetAssayData(visium_obj, layer="counts")[unique(random_gene_pairs %>% unlist()), ]

# For each gene pair,
null_dist <- lapply(1:n_random, function(i) {
  gene1 <- random_gene_pairs$gene1[i]
  gene2 <- random_gene_pairs$gene2[i]
  
  ligand_points <- tissue_coords %>%
    .[summary(random_genes_mat[gene1, ,drop=FALSE])$j, ] %>% 
    st_as_sf(coords = c("x", "y"))
  
  receptor_points <- tissue_coords %>%
    .[summary(random_genes_mat[gene2, ,drop=FALSE])$j, ] %>%
    st_as_sf(coords = c("x", "y"))
  
  # For each sender, find the nearest receiver
  lig_to_rec <- st_nearest_feature(ligand_points, receptor_points)
  
  # Calculate distance to nearest receptor
  ligand_points$dist_to_receptor <- st_distance(ligand_points, receptor_points[lig_to_rec, ], by_element = TRUE)
  
  ligand_points %>%
    mutate(dist_to_receptor_um = as.numeric(dist_to_receptor) * (bin_size / square_size)) %>% 
    pull(dist_to_receptor_um)
})

# saveRDS(null_dist, file = paste0("visium_hd_liver", dataset, "null_dist_1000gene_distance_", bin_size_str, ".rds"))

null_dist <- readRDS(file = paste0("visium_hd_liver", dataset, "null_dist_1000gene_distance_", bin_size_str, ".rds"))

null_dist_df <- null_dist %>% reshape2::melt() %>% 
  `colnames<-`(c("distance_um", "pair")) %>% 
  mutate(distance_um_round = round(distance_um))

null_dist_df %>%
  filter(distance_um_round < 500) %>%
  ggplot(aes(x = distance_um_round, group=pair)) +
  geom_histogram(binwidth = 8, alpha=0.1, position="identity", color="grey")

# TODO: DO THIS WITH RANDOM LIGAND-RECEPTOR PAIRS, INSTEAD OF RANDOM GENE PAIRS
# Read NicheNet ligand-receptor pairs
lr_pairs <- readRDS("~/Documents/nichenet/nichenet_v2/lr_network_mouse_21122021.rds")

lr_pairs <- lr_pairs %>%
  filter(from %in% rownames(visium_obj) & to %in% rownames(visium_obj))

n_random <- 1000
random_gene_pairs <- tibble(
  ligand = sample(lr_pairs$from, n_random, replace = TRUE),
  receptor = sample(lr_pairs$to, n_random, replace = TRUE)
)

all(!random_gene_pairs$ligand == random_gene_pairs$receptor)

# Get all info 
tissue_coords <- GetTissueCoordinates(visium_obj)
random_genes_mat <- GetAssayData(visium_obj, layer="counts")[unique(random_gene_pairs %>% unlist()), ]

# For each gene pair
null_dist_lr <- lapply(1:n_random, function(i) {
  print(i)
  gene1 <- random_gene_pairs$ligand[i]
  gene2 <- random_gene_pairs$receptor[i]
  
  ligand_points <- tissue_coords %>%
    .[summary(random_genes_mat[gene1, ,drop=FALSE])$j, ] %>% 
    st_as_sf(coords = c("x", "y"))
  
  receptor_points <- tissue_coords %>%
    .[summary(random_genes_mat[gene2, ,drop=FALSE])$j, ] %>%
    st_as_sf(coords = c("x", "y"))
  
  # For each sender, find the nearest receiver
  lig_to_rec <- st_nearest_feature(ligand_points, receptor_points)
  
  # Calculate distance to nearest receptor
  ligand_points$dist_to_receptor <- st_distance(ligand_points, receptor_points[lig_to_rec, ], by_element = TRUE)
  
  ligand_points %>%
    mutate(dist_to_receptor_um = as.numeric(dist_to_receptor) * (bin_size / square_size)) %>% 
    pull(dist_to_receptor_um)
})

saveRDS(null_dist_lr, file = paste0("visium_hd_liver", dataset, "null_dist_1000LR_distance_", bin_size_str, ".rds"))

null_dist_lr <- readRDS(file = paste0("~/visium-hd/visium_hd_liver", dataset, "null_dist_1000LR_distance_", bin_size_str, ".rds"))
square_size <- visium_obj@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot

null_dist_lr_df <- null_dist_lr %>% reshape2::melt() %>% 
  `colnames<-`(c("distance_um", "pair")) %>% 
  mutate(distance_um_round = round(distance_um))

null_dist_lr_df %>%
  # filter(distance_um_round < 500) %>%
  ggplot(aes(x = distance_um_round, group=pair)) +
  geom_histogram(binwidth = 8, alpha=0.1, position="identity", color="grey")

# Calculate distances between Gdf2 and Acvrl1
# Random gene 1
random_gene_2 <- sample(rownames(visium_obj), 2)


# Randomly shuffle Gdf2 labels
sample(GetAssayData(visium_obj, layer="counts")["Gdf2", ]) %>% names %>% head
GetAssayData(visium_obj, layer="counts")["Gdf2", ] %>% names %>% head

ligand_points_df <- lapply(c("Gdf2", "Bmp10"), function(gene){
  # ligand_points <- GetTissueCoordinates(visium_obj) %>%
  #   .[summary(genes_oi_mat[gene, ,drop=FALSE])$j, ] %>% 
  #   st_as_sf(coords = c("x", "y"))
  gene <- random_gene_2[1]
  gene_oi_mat <- GetAssayData(visium_obj, layer="counts")[random_gene_2, ]
  
  all_points <- GetTissueCoordinates(visium_obj) %>% 
    st_as_sf(coords = c("x", "y"))
  
  receptor_points <- GetTissueCoordinates(visium_obj) %>%
    .[summary(gene_oi_mat[random_gene_2[2] ,,drop=FALSE])$j, ] %>%
    st_as_sf(coords = c("x", "y"))
  
  # For each sender, find the nearest receiver
  lig_to_rec <- st_nearest_feature(all_points, receptor_points)
  
  # Calculate distance to nearest receptor
  all_points$dist_to_receptor <- st_distance(all_points, receptor_points[lig_to_rec, ], by_element = TRUE)
  
  tmp <- all_points %>%
    mutate(dist_to_receptor_um = as.numeric(dist_to_receptor) * (bin_size / square_size)) %>% 
    # join with count data
    bind_cols(gene_oi_mat[gene, .$cell] %>% data.frame(count=.)) %>% 
    mutate(gene = gene)

}) %>% bind_rows()


# Some tests
car::qqPlot(tmp %>% filter(count == 0) %>% pull(dist_to_receptor_um))
tmp %>% filter(count == 0 ) %>%  sample_n(5000) %>%  pull(dist_to_receptor_um) %>% shapiro.test()
tmp1 <- tmp %>% filter(count == 0) %>% pull(dist_to_receptor_um)
tmp2 <- tmp %>% filter(count > 0) %>% pull(dist_to_receptor_um)
wilcox.test(tmp1, tmp2)
t.test(tmp1, tmp2)

null_gene_label_shuffling <- sapply(1:1000, function(i){
  tmp_shuffled <- tmp %>%
    mutate(count = sample(count))
  
  mean_distance_expressed <- tmp_shuffled %>%
    filter(count > 0) %>%
    pull(dist_to_receptor_um) %>%
    mean
  
  return(mean_distance_expressed)
})

tmp2 <- tmp %>% filter(count > 0)

1 - (sum(mean(tmp2$dist_to_receptor_um) <= null_gene_label_shuffling) / length(null_gene_label_shuffling))
null_gene_label_shuffling

tmp %>%
  mutate(expressed = count > 0) %>%
  group_by(expressed) %>% 
  summarise(mean_distance = mean(dist_to_receptor_um),
            median_distance = median(dist_to_receptor_um),
            n = n())



ggplot(tmp_random, aes(x = dist_to_receptor_um, y=count)) +
  geom_jitter(size=0.5, stroke=0, shape=16) +
  facet_wrap(~gene, scales = "free") +
  labs(y = paste0("Raw counts in bins expressing the ligand"),
       x = paste0("Distance to nearest bin expressing Acvrl1 (\u03bcm)")) +
  theme_classic(base_size=7) +
  theme(axis.title.y = element_text(margin = margin(r=6))
  )

ggsave(paste0("plots/liver_ligand_receptor_distance_vs_counts_", bin_size_str, ".pdf"),
       device = cairo_pdf,
       width = 5, height = 3)

mean(ligand_points_df$dist_to_receptor_um)

# How many points are neighbors?
ligand_points_df %>% 
  filter(dist_to_receptor_um <= 12) %>% 
  group_by(gene) %>%
  summarise(n = n())

# Calculate mean distance for each permutation in the null dist
null_dist_lr_means <- sapply(null_dist_lr, mean, na.rm=TRUE) %>% 
  # omit NaNs
  .[!is.nan(.)]
# Number of times the mean distance in the null is less than or equal to the observed mean distance
p_value <- sum(null_dist_lr_means <= mean(ligand_points_df$dist_to_receptor_um)) / length(null_dist_lr_means)

# p_value of Gdf2-Acvrl1 = 0.406438
# p_value of Bmp10-Acvrl1 = 0.40744 

ggplot(ligand_points_df, aes(x = dist_to_receptor_um)) +
  geom_histogram(binwidth = 8)
