# Convert h5ad to seurat
# I had to do this in the terminal
# conda activate python_r, then launch R
# library(sceasy)
# library(reticulate)
# use_condaenv('python_r')
# 
# file_name <- "~/Documents/other_scripts/vizgen_data"
# sceasy::convertFormat(paste0(file_name, ".h5ad"),
#                       from="anndata", to="seurat",
#                       main_layer = "data",
#                       outFile=paste0(file_name, ".rds"))

library(Seurat)
library(tidyverse)
library(precrec)

seuratObj <- readRDS("~/Documents/other_scripts/vizgen_data.rds")

# Get cluster labels
metadata <- read.csv("~/Documents/other_scripts/vizgen_cluster_labels.csv", row.names = 1) %>% 
  cbind(., seuratObj@reductions$spatial@cell.embeddings) %>% 
  cbind(., seuratObj@reductions$umap@cell.embeddings) %>% 
  mutate(cell_id = 1:n())

seuratObj <- AddMetaData(seuratObj, metadata = metadata)

ggplot(metadata, aes(x=UMAP_1, y=UMAP_2, color=factor(leiden_covet_kcs))) +
  geom_point() +
  theme_bw()

ggplot(metadata %>% arrange(!is.na(leiden_covet_kcs)),
       aes(x=SPATIAL_1, y=SPATIAL_2, color=factor(leiden_covet_kcs))) +
  geom_point() + 
  scale_color_discrete(na.value = "gray90") +
  theme_bw()

ggplot(metadata %>% arrange(!is.na(neighborhood_cluster_kcs)),
       aes(x=SPATIAL_1, y=SPATIAL_2, color=factor(neighborhood_cluster_kcs))) +
  geom_point() + 
  scale_color_discrete(na.value = "gray90") +
  theme_bw()

# No zonation genes in the object
c('Glul', 'Cyp2e1') %in% rownames(seuratObj)
c('Ass1', 'Asl10', 'Alb8', 'Cyp2f2') %in% rownames(seuratObj) 

# How many ligands are in the data?
lr_network <- readRDS("~/Documents/nichenet/nichenet_v2/lr_network_mouse_21122021.rds")
ligand_target_matrix <- readRDS("~/Documents/nichenet/nichenet_v2/ligand_target_matrix_nsga2r_final_mouse.rds")

ligand_target_matrix %>% reshape2::melt()

ligand_target_matrix %>% reshape2::melt() %>% rename(target=Var1, ligand=Var2) %>% 
  ggplot(aes(x=value, group=ligand)) +
  geom_density(alpha=0.01) +
  theme_classic()

# Check overlap - vizgen vs resolve
intersect(lr_network %>% pull(from)  %>% unique, rownames(seuratObj)) %>% length # 59 ligands
intersect(lr_network %>% pull(to) %>% unique, rownames(seuratObj)) %>% length # 70 receptors
intersect(lr_network %>% select(from, to) %>% unlist %>% unique, rownames(seuratObj)) %>% length # 112 together

ligands_list <- intersect(lr_network %>% pull(from)  %>% unique, rownames(seuratObj))

# Get 8 nearest neighbors of each cell
seuratObj <- FindNeighbors(seuratObj, reduction="spatial", k.param=8, dims=1:2)

# Or do this manually to be able to set max distance
nns <- RANN::nn2(seuratObj@reductions$spatial@cell.embeddings, k=9)
# Plot distance distribution of 8 nearest neighbors
nns$nn.dists[,-1] %>% rowMeans %>% data.frame(val=.) %>% ggplot(aes(x=val)) + geom_density()
q_radius <- nns$nn.dists[,-1] %>% rowMeans %>% quantile(0.95)
# Set cutoff at 95th quantile
nns <-  RANN::nn2(seuratObj@reductions$spatial@cell.embeddings, k=9,
                  searchtype='radius', radius = q_radius)

# Plot ligand expression in mean of microenvironment vs distribution across all cells in seurat obj
plots <- lapply(ligands_list, function(lig) {
  ligand_expr <- as.vector(seuratObj[lig, ]@assays$RNA@data)
  nn_idx_df <- nns$nn.idx %>% data.frame() %>% pivot_longer(X2:X9) %>%
    select(-name) %>% rename(cell_id = X1, neighbor_id = value) %>%
    filter(neighbor_id != 0) %>%
    mutate(ligand = ligand_expr[neighbor_id])

  rbind(nn_idx_df %>% group_by(cell_id) %>%
          summarise(mean_lig = max(ligand)) %>%
          mutate(group='neighbor'),
        ligand_expr %>% data.frame(mean_lig = .) %>%
          mutate(cell_id = 1:n(), group='expr')) %>%
    ggplot(aes(x=mean_lig, color=group)) + geom_density() +
    theme_classic() + labs(title=lig)
})

names(plots) <- ligands_list

ggsave(filename = "plots/vizgen/max_expr_vs_null_1.png",
       plot = patchwork::wrap_plots(plots[1:30], ncol=6, nrow=5) &
  theme(legend.position = "none",
        axis.title = element_blank()),
  width = 20, height=12)

ggsave(filename = "plots/vizgen/max_expr_vs_null_2.png",
       patchwork::wrap_plots(plots[31:59], ncol=6, nrow=5) &
  theme(legend.position = "none",
        axis.title = element_blank(),
        width=20, height=12))

lig <- ligands_list[1]
ligand_expr <- as.vector(seuratObj[lig, ]@assays$RNA@data)
#%>% 
  #mutate(ligand = ligand_expr[neighbor_id])

# top 250 of lt_matrix
q_lt <- ligand_target_matrix[, lig] %>% sort(decreasing = TRUE) %>% .[1:250] %>% names
lt_vector <- ligand_target_matrix[rownames(ligand_target_matrix) %in% rownames(seuratObj), lig] %>% 
  replace(., 1:length(.), 0) %>% .[sort(names(.))]
lt_vector[names(lt_vector) %in% q_lt] <- 1

obs <- seuratObj[rownames(seuratObj) %in% rownames(ligand_target_matrix), 1]@assays$RNA@data

c(obs)

prcs <- precrec::evalmod(scores=c(obs), labels=lt_vector, mode="prcroc")
attr(prcs$prcs[[1]], "auc")

#### NICHENET-LIKE ####
# Get a dataframe of a cell and all its neighbors
nn_idx_df <- nns$nn.idx %>% data.frame() %>% pivot_longer(X2:X9) %>% 
  select(-name) %>% rename(cell_id = X1, neighbor_id = value) %>% 
  filter(neighbor_id != 0) 

cell_ids_list <- nn_idx_df$cell_id %>% unique %>% .[1:10000] # Subset

# Subset objects for faster lookup
seuratObj_subset <- seuratObj[rownames(seuratObj) %in% rownames(ligand_target_matrix),]@assays$RNA@data
seurat_expr_ligands <- seuratObj_subset[ligands_list,]
lt_matrix_subset <- ligand_target_matrix[rownames(ligand_target_matrix) %in% rownames(seuratObj), ]

# For each cell
prcs_list <- pbmcapply::pbmclapply(cell_ids_list, function(id) {
  #print(id)
  
  # Get nearest neighbors of each cell
  neighbor_ids <- nn_idx_df %>% filter(cell_id == id) %>% pull(neighbor_id)
  
  # Get expressed ligands from microenvironment (if any of the neighbors express the ligand)
  expressed_ligands <- rowSums(seurat_expr_ligands[, neighbor_ids, drop=FALSE] > 0) %>%
    .[. > 0] %>% names
  
  # TODO: FILTER BASED ON RECEPTOR EXPRESSION IN RECEIVER CELL ?
  
  # Get expression of the receiver cell
  obs_expr <- seuratObj_subset[, id]
  
  # For each expressed ligand, create a ground truth vector of target genes
  # 1 = gene is in the top 250 target of the ligand (over the whole L-T matrix)
  labels <- lapply(expressed_ligands, function(expr_lig) {
    # Get top 250 of the ligand in the LT matrix
    q_lt <- ligand_target_matrix[, expr_lig] %>% sort(decreasing = TRUE) %>% .[1:250] %>% names
    
    # Create vector of zeroes -> sort alphabetically
    lt_vector <- lt_matrix_subset[, expr_lig] %>% 
      replace(., 1:length(.), 0) %>% .[sort(names(.))]
    
    # Set to 1 if gene is in the top 250
    lt_vector[names(lt_vector) %in% q_lt] <- 1
    
    lt_vector %>% unname
  })
  
  # Create score matrix based on real expression
  scores <- rep(list(c(obs_expr)), length(labels))
  
  # Calculate AUPR for each ligand
  prcs <- precrec::evalmod(scores=scores, labels=labels,
                           dsids = 1:length(expressed_ligands), mode="prcroc")
  prcs_cell <- attr(prcs, "aucs") %>% filter(curvetypes == "PRC")
  prcs_cell$aucs %>% setNames(expressed_ligands)
  
}, mc.cores = 12)

# There was some error in prcs_list -> do this for now
prcs_list_c <- prcs_list[which(sapply(prcs_list, is.double))] %>%
  setNames(which(sapply(prcs_list, is.double)))

# Plot distribution of max AUPR
sapply(prcs_list_c, max) %>% data.frame(val=.) %>% ggplot(aes(x=val)) + geom_density() + theme_classic()

# Add AUPR to metadata, plot distribution spatially
metadata$aupr <- NA
metadata[as.numeric(names(prcs_list_c)), "aupr"] <- sapply(prcs_list_c, max)
ggplot(metadata %>% arrange(!is.na(aupr)),
       aes(x=SPATIAL_1, y=SPATIAL_2, color=aupr)) +
  geom_point() + 
  scale_color_viridis_b(na.value = "gray90") +
  theme_bw()

# Plotting a specific ligand expression
loi <- "Angpt2"
FeaturePlot(seuratObj, loi, reduction="spatial", order=TRUE)

# Calculate Moran's I
dists_matrix <- as.matrix(dist(metadata[as.numeric(names(prcs_list_c)), c("SPATIAL_1", "SPATIAL_2")]))
dists_matrix <- 1/dists_matrix
diag(dists_matrix) <- 0
moran <- ape::Moran.I(sapply(prcs_list_c, max), dists_matrix)
moran

# Investigating results
# Melt into dataframe, keep names of ligands
prcs_df <- prcs_list_c %>% unlist %>% data.frame(aupr=.) %>% 
  rownames_to_column("cell_ligand") %>% 
  separate(cell_ligand, c("cell_id", "ligand"), sep="\\.")

# Join results with cell type
prcs_df_celltype <- left_join(prcs_df %>% mutate(cell_id = as.numeric(cell_id)),
                              seuratObj$annotation_own_score_genes %>% data.frame(celltype=.) %>%
                                mutate(cell_id = 1:n()),
                              by="cell_id") 

# Summarise prcs_df_celltype_sum per celltype
prcs_df_celltype_summ <- prcs_df_celltype %>%  group_by(celltype, ligand) %>%
  summarise(mean_aupr = median(aupr)) %>% arrange(desc(mean_aupr))

prcs_df_celltype_summ


### CODE DUMP ####
# For each ligand, get the expression in the 8 nearest neighbors
lapply(ligands_list, function(ligand) {
  # Replace ones in seuratObj@graphs$RNA_nn with the actual gene expression
  
})

Idents(seuratObj) <- seuratObj$leiden_covet
Idents(seuratObj) <- ifelse(seuratObj$leiden_covet == 1, "1", "other")

markers <- FindMarkers(seuratObj,
               features = ligands_list,
               ident.1 = "1",
               ident.2 = "other",
               min.pct = 0.05,
               only.pos = TRUE,
               return.thresh = 1,
               logfc.threshold = 0)

which(seuratObj@graphs$RNA_nn == 1)


# dgTMatrix is a lot easier to work with, but still smaller than a dense matrix
neighbor_graph_T_default <- as(neighbor_graph, "TsparseMatrix")

# get_neighbor_expr <- sapply(0:max(neighbor_graph_T@i), function(i){
#   indices <- neighbor_graph_T@j[which(neighbor_graph_T@i == i)] + 1
#   as(seuratObj[ligand, indices]@assays$RNA@counts, "TsparseMatrix")
#   
# })

neighbor_graph_T <- neighbor_graph_T_default
ligand <- "Nrp1"
neighbor_graph_T@x <- rep(as.vector(seuratObj[ligand,]@assays$RNA@data), times= table(neighbor_graph_T@j))


neighbor_df <- data.frame(row=neighbor_graph_T@i, col=neighbor_graph_T@j, val=neighbor_graph_T@x) %>% 
  arrange(row) %>% mutate(cell_id = row + 1) %>% 
  # Join metadata
  left_join(metadata %>% mutate(cell_id = 1:n()) %>%
              select(cell_id, leiden_covet, annotation_own_score_genes),
            by="cell_id")

ggplot(neighbor_df, aes(x=factor(leiden_covet), y=val)) +
  geom_boxplot() +
  theme_bw()

sapply(0:9, function(i){
  indices <- neighbor_graph_T@j[which(neighbor_graph_T@i == i)] + 1
  as(seuratObj[ligand, indices]@assays$RNA@counts, "TsparseMatrix")
})

neighbor_transposed <- Matrix::t(neighbor_graph_T)
neighbor_new <- Matrix::drop0(neighbor_transposed, give.Csparse = FALSE)

which(neighbor_graph[1, ] == 1)