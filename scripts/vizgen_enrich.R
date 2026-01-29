library(tidyverse)
library(Giotto)
library(patchwork)
first_run <- FALSE

if (first_run){
  metadata <- read.csv("data/vizgen_annotation.csv", row.names = 1) %>% 
    rename(celltype = annotation_own_score_genes)
  expr <- read.csv("data/vizgen_polyT_counts.csv")
  expr <- expr %>% column_to_rownames("cells") %>% t
  all(rownames(metadata) == colnames(expr))
  
  giotto_obj <- createGiottoObject(expr, spatial_locs = metadata %>% select(x, y),
                                         cell_metadata = metadata %>% select(celltype) %>% 
                                     rownames_to_column("cell_ID"))
  
  giotto_obj <- createSpatialDelaunayNetwork(gobject = giotto_obj)
  giotto_obj <- createSpatialKNNnetwork(gobject = giotto_obj, k=5)
  
  saveRDS(giotto_obj, "data/vizgen/vizgen_giotto_obj.rds")
} else {
  giotto_obj <- readRDS("data/vizgen/vizgen_giotto_obj.rds")
}

# Get roi 
roi <- c(10000,15000,8000,12000)


giotto_obj_subset <- subsetGiottoLocs(gobject = giotto_obj,
                                      x_min = roi[1],
                                      x_max = roi[2],
                                      y_min = roi[3],
                                      y_max = roi[4])
spatPlot(giotto_obj_subset) + scale_y_reverse()
spatPlot(gobject = giotto_obj_subset,
         show_network = T,         point_shape = "no_border",
         network_color = 'black',
         spatial_network_name = 'Delaunay_network',
         point_size = 1.5,
         cell_color = "celltype",
         show_plot = FALSE,
         coord_fix_ratio = 1) + scale_y_reverse()

spatPlot(gobject = giotto_obj_subset,
         show_network = T,
         point_shape = "no_border",
         network_color = 'black',
         spatial_network_name = 'knn_network',
         point_size = 1.5,
         cell_color = "celltype",
         coord_fix_ratio = 1) + scale_y_reverse()

plotStatDelaunayNetwork(gobject = giotto_obj)

# Create other spatial networks
giotto_obj_maxdistNULL <- createSpatialDelaunayNetwork(gobject = giotto_obj, maximum_distance = NULL)
saveRDS(giotto_obj_maxdistNULL@spatial_network$cell$Delaunay_network, "rds/vizgen_giotto_delaunay_network_maxdistNULL.rds")

for (k in c(5, 7, 10)) {
  giotto_obj <- createSpatialKNNnetwork(gobject = giotto_obj, k=k)
  giotto_obj_subset <- subsetGiottoLocs(gobject = giotto_obj,
                                        x_min = roi[1],
                                        x_max = roi[2],
                                        y_min = roi[3],
                                        y_max = roi[4])
  
  p_knn <- spatPlot(gobject = giotto_obj_subset,
                    show_network = T,
                    point_shape = "border",
                    network_color = 'gray70',
                    point_border_col = "black",
                    point_border_stroke = 0.2,
                    spatial_network_name = 'knn_network',
                    point_size = 1.5,
                    cell_color = "celltype",
                    cell_color_code = color_palette,
                    show_plot = FALSE,
                    coord_fix_ratio = 1) +
    scale_y_reverse() +
    ggtitle (paste0("KNN graph with k=",k)) +
    theme_void() +
    theme(panel.background = element_rect(fill = "white"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 8, face = "bold"))
  
  # Change linewidth of the segment
  p_knn$layers[[1]]$aes_params$linewidth <- 0.2
  
  ggsave(paste0("plots/vizgen/vizgen_knn_network_k", k, ".png"),
         plot = p_knn,
         width = 4,height = 4,
         dpi = 300)
}

# Calculate cell proximities
spatial_network <- 'knn_network'
for (spatial_network in c('Delaunay_network', 'knn_network')) {
  
  if (spatial_network == 'knn_network') {
    
    for (k in c(5,7,10)){
      giotto_obj <- createSpatialKNNnetwork(gobject = giotto_obj, k=5)
      
      cell_proximities <- cellProximityEnrichment(gobject = giotto_obj,
                                                  cluster_column = 'celltype',
                                                  spatial_network_name = spatial_network,
                                                  adjust_method = 'fdr',
                                                  number_of_simulations = 1000)
      saveRDS(cell_proximities, paste0("rds/vizgen_giotto_cell_proximities_knn_network_k",
                                         k, ".rds"))
    }
                                       
  }
  else {
    cell_proximities <- cellProximityEnrichment(gobject = giotto_obj,
                                                cluster_column = 'celltype',
                                                spatial_network_name = spatial_network,
                                                adjust_method = 'fdr',
                                                number_of_simulations = 1000)
    
    saveRDS(cell_proximities, paste0("rds/vizgen_giotto_cell_proximities_delaunay.rds"))
    
  }
}


cell_proximities <- readRDS("vizgen_giotto_cell_proximities_delaunay.rds")


cell_proximities_list <- list(
  delaunay = readRDS("rds/vizgen_giotto_cell_proximities_delaunay.rds"),
  knn_k5 = readRDS("rds/vizgen_giotto_cell_proximities_knn_k5.rds"),
  knn_k7 = readRDS("rds/vizgen_giotto_cell_proximities_knn_k7.rds"),
  knn_k10 = readRDS("rds/vizgen_giotto_cell_proximities_knn_k10.rds")
)

cell_proximities

cellProximityBarplot(gobject = giotto_obj,
                     CPscore = cell_proximities_list[[3]],
                     min_orig_ints = 3,
                     min_sim_ints = 3)

cellProximityNetwork(gobject = giotto_obj,
                     CPscore = cell_proximities,
                     remove_self_edges = F,
                     self_loop_strength = 0.3,
                     only_show_enrichment_edges = T,
                     rescale_edge_weights = T,
                     node_size = 8,
                     edge_weight_range_depletion = c(1,2),
                     edge_weight_range_enrichment = c(2,5))

p_heatmaps <- lapply(cell_proximities_list, 
                     function(i) {cellProximityHeatmap(gobject = giotto_obj,
                     CPscore = i,
                     order_cell_types = FALSE,
                     scale=FALSE)})


cellProximityHeatmapMatrix <- function(
    gobject,
    CPscore) {
  enrich_res <- CPscore$enrichm_res
  
  # data.table variables
  first_type <- second_type <- unified_int <- NULL
  
  enrich_res[, first_type := strsplit(
    x = as.character(unified_int), split = "--"
  )[[1]][1],
  by = seq_len(nrow(enrich_res))
  ]
  enrich_res[, second_type := strsplit(
    x = as.character(unified_int), split = "--"
  )[[1]][2],
  by = seq_len(nrow(enrich_res))
  ]
  
  # create matrix
  enrich_mat <- data.table::dcast.data.table(
    data = enrich_res,
    formula = first_type ~ second_type,
    value.var = "enrichm"
  )
  matrix_d <- as.matrix(enrich_mat[, -1])
  rownames(matrix_d) <- as.vector(enrich_mat[[1]])
  t_matrix_d <- t_flex(matrix_d)
  
  # fill in NAs based on values in upper and lower matrix triangle
  t_matrix_d[upper.tri(t_matrix_d)][is.na(t_matrix_d[
    upper.tri(t_matrix_d)
  ])] <- matrix_d[upper.tri(matrix_d)][
    is.na(t_matrix_d[upper.tri(t_matrix_d)])
  ]
  t_matrix_d[lower.tri(t_matrix_d)][is.na(t_matrix_d[
    lower.tri(t_matrix_d)
  ])] <- matrix_d[lower.tri(matrix_d)][
    is.na(t_matrix_d[lower.tri(t_matrix_d)])
  ]
  t_matrix_d[is.na(t_matrix_d)] <- 0
  final_matrix <- t_matrix_d
  
  return (final_matrix)
}

for(name in names(cell_proximities_list)) {
  heatmap_matrix <- cellProximityHeatmapMatrix(giotto_obj, cell_proximities_list[[name]])

  saveRDS(heatmap_matrix,
          paste0("rds/vizgen_giotto_cell_proximities_heatmap_matrix_", name, ".rds"))
}


cell_proximities$enrichm_res %>% filter(grepl("Kuppfer cells", unified_int))
cell_proximities$enrichm_res %>% filter(grepl("LAM", unified_int))

