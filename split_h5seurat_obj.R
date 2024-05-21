library(hdf5r)
library(Matrix)
library(tidyverse)

library(Seurat)

# Read in h5Seurat object
seurat_obj <- H5File$new("LiverDataReleaseSeurat_newUMAP.h5Seurat", mode = "r")

# Explore this object a little bit
seurat_obj
unique(seurat_obj[["meta.data"]][["cellType"]][])
unique(seurat_obj[["meta.data"]][["Run_Tissue_name"]][]) # "NormalLiver" and "CancerousLiver"
unique(seurat_obj[["meta.data"]][["fov"]][])
table(seurat_obj[["meta.data"]][["fov"]][], seurat_obj[["meta.data"]][["Run_Tissue_name"]][]) # Cancer has more FOVs

# Get the count matrix of the object
hdf5_counts <- seurat_obj[["assays"]][["RNA"]][["counts"]]
countmatrix <- sparseMatrix(i = hdf5_counts[["indices"]][],
                     x = hdf5_counts[["data"]][],
                     p = hdf5_counts[["indptr"]][],
                     index1 = FALSE)
colnames(countmatrix) <- seurat_obj[["cell.names"]][]
rownames(countmatrix) <- seurat_obj[["assays"]][["RNA"]][["features"]][]

# Get the normalized count matrix of the object
hdf5_counts_norm <- seurat_obj[["assays"]][["RNA"]][["data"]]
norm_countmatrix <- sparseMatrix(i = hdf5_counts_norm[["indices"]][],
                            x = hdf5_counts_norm[["data"]][],
                            p = hdf5_counts_norm[["indptr"]][],
                            index1 = FALSE)
colnames(norm_countmatrix) <- seurat_obj[["cell.names"]][]
rownames(norm_countmatrix) <- seurat_obj[["assays"]][["RNA"]][["features"]][]

# Dimensionality reductions
names(seurat_obj[["reductions"]])

# Only get UMAP
reduc_df <- lapply(names(seurat_obj[["reductions"]])[1:2], function(reduc_col){
  print(reduc_col)
  seurat_obj[["reductions"]][[reduc_col]][["cell.embeddings"]][,] %>%
    `colnames<-`(paste0(reduc_col, "_", 1:2))
}) %>% do.call(cbind, .) %>% data.frame %>% mutate(cell_ID = seurat_obj[["cell.names"]][])


# Check metadata column names
h5attr(seurat_obj[["meta.data"]], "colnames")
metadata_cols <- h5attr(seurat_obj[["meta.data"]], "colnames")
metadata_cols

cols_convert <- c("fov", "x_FOV_px", "y_FOV_px", "x_slide_mm", "y_slide_mm", "idx") # convert to numeric
metadata_df <- sapply(metadata_cols, function(col){
  seurat_obj[["meta.data"]][[col]][]
}) %>% data.frame() %>% rownames_to_column("idx") %>%
  mutate(across(all_of(cols_convert), ~as.numeric(.x)))

# Join metadat and UMAP information
metadata_df <- metadata_df %>% inner_join(reduc_df, by = "cell_ID")

# Create seurat object
seurat_obj_all <-  CreateSeuratObject(counts = countmatrix,
                                      meta.data = metadata_df %>% column_to_rownames("cell_ID"))
saveRDS(seurat_obj_all, paste0("LiverData/LiverData_Raw.rds"))
# Add normalized data
seurat_obj_all <- SetAssayData(seurat_obj_all, slot = "data", new.data = norm_countmatrix)
saveRDS(seurat_obj_all, paste0("LiverData/LiverData_RawNorm.rds"))

# Check indices of normal and cancer
R.utils::seqToIntervals(grep("Normal", seurat_obj_all$Run_Tissue_name))
R.utils::seqToIntervals(grep("Cancerous", seurat_obj_all$Run_Tissue_name))

# Subset each run into four quadrants
for (run_oi in c("NormalLiver", "CancerousLiver")){
  seurat_obj_subset <- seurat_obj_all[, seurat_obj_all$Run_Tissue_name == run_oi]
  gc()
  meta_df_subset <- seurat_obj_subset@meta.data
  
  # Get mean position for each FOV
  text_pos <- meta_df_subset %>% group_by(fov) %>% summarise(mean_x = mean(x_slide_mm), mean_y = mean(y_slide_mm))
  
  # Get quadrants
  x_mid <- mean(c(max(meta_df_subset$x_slide_mm), min(meta_df_subset$x_slide_mm)))
  y_mid <- mean(c(max(meta_df_subset$y_slide_mm), min(meta_df_subset$y_slide_mm)))
  
  # Label each cell per quadrant
  meta_df_subset <- meta_df_subset %>% 
    mutate(quadrant = case_when(x_slide_mm > x_mid & y_slide_mm > y_mid   ~ "Q1",
                                x_slide_mm <= x_mid & y_slide_mm > y_mid  ~ "Q2",
                                x_slide_mm <= x_mid & y_slide_mm <= y_mid ~ "Q3",
                                TRUE                                      ~ "Q4"))
  
  seurat_obj_subset$quadrant <- meta_df_subset$quadrant
  
  p <- ggplot(seurat_obj_subset@meta.data, aes(x = x_slide_mm, y = y_slide_mm, color = quadrant)) +
    geom_point(size=0.2) +
    geom_text(data=text_pos, aes(label = fov, x=mean_x, y = mean_y), inherit.aes = FALSE) +
    theme_classic() +
    geom_vline(xintercept = x_mid) + # plot vertical line
    geom_hline(yintercept = y_mid) + coord_fixed()
  ggsave(paste0("LiverData/quadrant_plot_", run_oi, ".png"), p,
         width = 3000, height = 3000, units = "px")
  
  # Subset into quadrants
  for (q_oi in paste0("Q", 1:4)){
    seurat_obj_subset_q <- seurat_obj_subset[, seurat_obj_subset$quadrant == q_oi]
    saveRDS(seurat_obj_subset_q, paste0("LiverData/LiverData_RawNorm_", run_oi, "_", q_oi, ".rds"))
  }
  

}

# conversion
library(sceasy)
library(reticulate)

use_condaenv("squidpy")

liver_seurat_obj <- readRDS("LiverData/LiverData_RawNorm.rds")
sceasy::convertFormat(liver_seurat_obj, from="seurat", to="anndata",
                      outFile='LiverData_RawNorm.h5ad')




