library(Seurat)
library(tidyverse)
library(Matrix)

data_path_root <- "~/visium-hd/data/"
data_paths <- list(
  liver_visium_sca = "Visium_HD_Liver/Visium_HD_Liver_008um.rds",
  liver_visium_caw = "Visium_HD_Liver_CAW009/Visium_HD_Liver_CAW009_008um.rds",
  brain_visium_ffpe = "Visium_HD_MouseBrain/Visium_HD_MouseBrain_008um.rds",
  brain_visium_ff = "Visium_HD_MouseBrain_FF/Visium_HD_MouseBrain_FF_008um.rds",
  lung_visium = "Visium_HD_Lung_ILA010/Visium_HD_Lung_ILA010_008um.rds",
  liver_sc = "scref_Liver_Guilliams/liver_mouseStSt_noEC.rds",
  brain_sc = "scref_MouseBrain_ABA/WMB-10Xv3_subset.rds",
  lung_sc = "scref_Lung_UBla/lung_combined_macs_DCs_seurat.rds"
)

# My savior: https://rpubs.com/will_townes/sparse-apply
listCols<-function(m){
  # Converts a sparse Matrix into a list of its columns (only non-zeroes entries)
  res<-split(m@x, findInterval(seq_len(nnzero(m)), m@p, left.open=TRUE))
  names(res)<-colnames(m)
  res
}

dataset <- "liver_visium_sca"
dataset <- "brain_visium_ffpe"
dataset <- "lung_sc"

counts_per_genes_df <- lapply(names(data_paths), function(dataset){
  cat("Processing dataset:", dataset, "\n")
  seurat_obj <- readRDS(paste0(data_path_root, data_paths[[dataset]]))
  active_assay_name <- DefaultAssay(seurat_obj)
  
  if (dataset == "brain_sc"){
    seurat_obj$nCount_RNA <- GetAssayData(seurat_obj, layer = "counts") %>% colSums()
    seurat_obj$nFeature_RNA <- (GetAssayData(seurat_obj, layer = "counts") > 0) %>% colSums()
  } 
  
  # Get spots with zero counts
  zero_counts <- which(seurat_obj@meta.data[,paste0("nCount_", active_assay_name)] == 0)
  
  if (length(zero_counts)){
    cat("Warning: Some spots/cells have 0 counts, subsetting...\n")
    seurat_obj <- seurat_obj[, -zero_counts]
    gc()
  }
  

  nCount <- seurat_obj@meta.data[, paste0("nCount_", active_assay_name)]
  nFeature <- seurat_obj@meta.data[, paste0("nFeature_", active_assay_name)]
  
  # Get counts per gene on average
  avg_gene <- nCount / nFeature
  mm <- GetAssayData(seurat_obj, layer = "counts")
  rm(seurat_obj); gc()
  
  mean_test <- pbmcapply::pbmclapply(listCols(mm), mean, mc.cores = 4) %>% unlist
  sd_vals <- pbmcapply::pbmclapply(listCols(mm), sd, mc.cores = 4) %>% unlist
  
  cat(paste0("Sanity check mean: ", all.equal(unname(mean_test), avg_gene), "\n"))
  
  # Split data set into tissue, technology, and sample
  dataset_parts <- strsplit(dataset, "_")[[1]]
  
  df <- data.frame(mean = avg_gene,
             sd = sd_vals) %>% 
    mutate(dataset_full = dataset,
           tissue = dataset_parts[1],
           tech = dataset_parts[2],
           sample = dataset_parts[3])
  
}) %>% bind_rows

saveRDS(counts_per_genes_df, file = paste0("rds/mean_sd_counts_per_gene_df.rds"))
counts_per_genes_df <- readRDS(file = paste0("rds/mean_sd_counts_per_gene_df.rds"))


# Plot histogram of visium
# Factor visium data
counts_per_genes_df <- counts_per_genes_df %>% 
  mutate(dataset_full = factor(dataset_full,
        levels = c("liver_visium_sca", "brain_visium_ffpe", "lung_visium",
                  "liver_visium_caw", "brain_visium_ff", 
                  "liver_sc", "brain_sc", "lung_sc")))

dataset_names <- c("Liver SCA002", "Brain FFPE", "Lung",
                   "Liver CAW009", "Brain Fresh Frozen", 
                   "Liver", "Brain", "Lung") %>% 
  setNames(levels(counts_per_genes_df$dataset_full))

counts_per_genes_df <- counts_per_genes_df %>% pivot_longer(cols = c("mean", "sd"),
                                    names_to = "statistic",
                                    values_to = "x")

set1_colors <- RColorBrewer::brewer.pal(8, "Set1")[1:2]

# Histogram of Visium HD
ggplot(counts_per_genes_df %>% filter(tech == "visium"), aes(x = x, fill=statistic)) +
  geom_histogram(color="white", bins=50, alpha=0.5) +
  facet_wrap(~dataset_full, scales = "free",
             labeller = labeller(dataset_full = dataset_names)) +
  scale_y_continuous(labels = scales::label_number(scale_cut=scales::cut_short_scale())) +
  theme_classic(base_size = 7) +
  theme(axis.text = element_text(size = 6),
        axis.title = element_blank(),
        strip.text = element_text(size = 6),
  )

ggplot(counts_per_genes_df %>% filter(tech == "visium", statistic == "mean"), aes(x=dataset_full, y = x)) +
  #geom_histogram(color="white", bins=50, alpha=0.5) +
  geom_boxplot() +
  # facet_wrap(~dataset_full, scales = "free",
  #            labeller = labeller(dataset_full = dataset_names)) +
  # scale_y_continuous(labels = scales::label_number(scale_cut=scales::cut_short_scale())) +
  theme_classic(base_size = 7) +
  theme(axis.text = element_text(size = 6),
        axis.title = element_blank(),
        strip.text = element_text(size = 6),
  )


# Histogram of scRNA-seq
ggplot(counts_per_genes_df %>% filter(tech == "sc"), aes(x = x, fill=statistic)) +
  geom_histogram(color="white", bins=50, alpha=0.5) +
  facet_wrap(~dataset_full, scales = "free",
             labeller = labeller(dataset_full = dataset_names)) +
  scale_y_continuous(labels = scales::label_number(scale_cut=scales::cut_short_scale())) +
  theme_classic(base_size = 7) +
  theme(axis.text = element_text(size = 6),
        axis.title = element_blank(),
        strip.text = element_text(size = 6),
  )

ggplot(counts_per_genes_df %>% filter(statistic == "mean"), aes(x=dataset_full, y = x)) +
  #geom_histogram(color="white", bins=50, alpha=0.5) +
  geom_boxplot() +
  # facet_wrap(~dataset_full, scales = "free",
  #            labeller = labeller(dataset_full = dataset_names)) +
  # scale_y_continuous(labels = scales::label_number(scale_cut=scales::cut_short_scale())) +
  theme_classic(base_size = 7) +
  theme(axis.text = element_text(size = 6),
        axis.title = element_blank(),
        strip.text = element_text(size = 6),
  )



counts_per_genes_df %>% 
  group_by(dataset_full, statistic) %>%
  summarise(mean_counts_per_gene = mean(x, na.rm=TRUE),
            median_counts_per_gene = median(x, na.rm=TRUE),
            n = n())

max_counts_df <- lapply(names(data_paths), function(dataset){
  seurat_obj <- readRDS(paste0(data_path_root, data_paths[[dataset]]))
  
  maxCounts <- sparseMatrixStats::colMaxs(GetAssayData(seurat_obj, layer = "counts"))
  
  dataset_parts <- strsplit(dataset, "_")[[1]]
  
  data.frame(x = maxCounts) %>% 
    mutate(dataset_full = dataset,
           tissue = dataset_parts[1],
           tech = dataset_parts[2],
           sample = dataset_parts[3])
}) %>% bind_rows


counts_per_genes_df %>% filter(tech == "sc", tissue == "brain")

# Compare counts of the two datasets
ggplot(counts_features_df, aes(x = x, fill=dataset)) +
  geom_histogram(color="white", alpha=0.5, bins=50, position="identity", size=0.1) +
  geom_text(aes(x = text_x, y = text_y, color=dataset,
                label = paste0("Mean: ", mean_x, "\nn = ", n)),
            data = means_df, vjust = 1, hjust=0, size=2,
            show.legend = FALSE
            ) +
  facet_wrap(type ~ bin_size, scales = "free",
             labeller = labeller( .multi_line=FALSE)) +
  scale_fill_manual(values = set1_colors, name = "Dataset") +
  scale_color_manual(values = set1_colors) +
  # Use scientific notation for y
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(labels = scales::comma) +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom",
        # Add right spacing to legend title
        legend.title = element_text(margin = margin(r = 15), size=8),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 6),
        axis.title = element_blank(),
        strip.text = element_text(size = 6),
        )
ggsave(paste0(plot_path, "histogram_counts_features.pdf"), width = 8, height = 5)


