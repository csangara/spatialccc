library(Giotto)
library(tidyverse)
library(patchwork)
# Plot Giotto and Squidpy neighborhood enrichment results
color_palette <- c("Hepatocytes" = "#B4B5B5FF",
                   "Endothelial cells" = "#dca754",
                   "Cholangiocytes" = "#C61B84FF",
                   "Stromal cells" = "#79151e",
                   "Kuppfer cells" = "#5DA6DBFF",
                   "Other immunecells" = "#893A86FF",
                   "B cells" = "#9C7EBAFF",
                   "Unknown" = "#191919")


giotto_obj <- readRDS("data/vizgen/vizgen_giotto_obj.rds")

giotto_obj@spatial_network$cell$Delaunay_network_squidpy <- readRDS("rds/vizgen_giotto_delaunay_network_maxdistNULL.rds")

giotto_obj <- setSpatialNetwork(
  gobject = giotto_obj,
  spat_unit = NULL,
  name = "Delaunay_network_squidpy",
  x = readRDS("rds/vizgen_giotto_delaunay_network_maxdistNULL.rds")
)

# Get roi 
roi <- c(10000,15000,8000,12000)
roi <- c(11000, 14000, 9000, 11500)

giotto_obj_subset <- subsetGiottoLocs(gobject = giotto_obj,
                                      x_min = roi[1],
                                      x_max = roi[2],
                                      y_min = roi[3],
                                      y_max = roi[4])

spatPlot(giotto_obj_subset) + scale_y_reverse()

p_spat_connect_titles <- c("Giotto Delaunay graph", "Squidpy Delaunay graph", "Giotto KNN graph (k=5)") %>% 
  setNames(c('Delaunay_network', 'Delaunay_network_squidpy', 'knn_network'))
# Spatial connectivity plots
p_spat_connect_plots <- lapply(c('Delaunay_network_squidpy', 'Delaunay_network', 'knn_network'), function(net) {
  p <- spatPlot(gobject = giotto_obj_subset,
           show_network = T,
           point_shape = "border",
           network_color = 'gray70',
           point_border_col = "black",
           point_border_stroke = 0.1,
           spatial_network_name = net,
           point_size = 1,
           cell_color = "celltype",
           cell_color_code = color_palette,
           show_plot = FALSE,
           coord_fix_ratio = 1) +
    scale_y_reverse() +
    ggtitle(p_spat_connect_titles[[net]]) +
    theme_void() +
    theme(#panel.background = element_rect(fill = "white"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 7, face = "bold"))
  
  # Change linewidth of the segment
  p$layers[[1]]$aes_params$linewidth <- 0.15
  
  p
})


wrap_plots(p_spat_connect_plots)


ggsave("plots/vizgen/vizgen_delaunay_network_giotto.png",
       plot = p_delaunay,
       width = 4,height = 4,
       dpi = 300)



giotto_obj$celltype




plotStatDelaunayNetwork(gobject = giotto_obj)


cell_proximities <- readRDS("vizgen_cell_proximities_delaunay.rds")
spatial_network <- 'knn_network'
for (spatial_network in c('Delaunay_network', 'knn_network')) {
  cell_proximities <- cellProximityEnrichment(gobject = giotto_obj,
                                              cluster_column = 'celltype',
                                              spatial_network_name = spatial_network,
                                              adjust_method = 'fdr',
                                              number_of_simulations = 1000)
  
  
  saveRDS(cell_proximities, paste0("rds/vizgen_cell_proximities_",
                                   tolower(str_replace(spatial_network, "_network", "")),
                                   ".rds"))
  
}

cell_proximities_list <- list(
  delaunay = readRDS("rds/vizgen_giotto_cell_proximities_delaunay.rds"),
  knn_k5 = readRDS("rds/vizgen_giotto_cell_proximities_knn_k5.rds"),
  knn_k7 = readRDS("rds/vizgen_giotto_cell_proximities_knn_k7.rds"),
  knn_k10 = readRDS("rds/vizgen_giotto_cell_proximities_knn_k10.rds")
)

cell_proximities

cellProximityBarplot(gobject = giotto_obj,
                     CPscore = cell_proximities_list[[1]],
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

celltypes <- sort(unique(giotto_obj$celltype)) %>% 
  # Relace "Kuppfer cells" with "Kupffer cells"
  str_replace("Kuppfer cells", "Kupffer cells") %>% 
  # "Other immunecells" with "Other immune cells"
  str_replace("Other immunecells", "Other immune cells")

cell_proximity_matrix <- list(
  delaunay = readRDS("rds/vizgen_giotto_cell_proximities_heatmap_matrix_delaunay.rds"),
  knn_k5 = readRDS("rds/vizgen_giotto_cell_proximities_heatmap_matrix_knn_k5.rds"),
  knn_k7 = readRDS("rds/vizgen_giotto_cell_proximities_heatmap_matrix_knn_k7.rds"),
  knn_k10 = readRDS("rds/vizgen_giotto_cell_proximities_heatmap_matrix_knn_k10.rds"),
  squidpy = read.table("rds/vizgen_squidpy_nhood_enrichment_zscore.txt", sep="\t") %>% as.matrix()
) %>% lapply(function(x) {x %>% `rownames<-`(celltypes) %>% `colnames<-`(celltypes)})


lims_list <- list(
  delaunay = seq(-4, 4, length = 3),
  knn_k5 = seq(-4, 4, length = 3),
  squidpy = seq(-100, 100, length = 3)
)

p_heatmaps <- lapply(c("squidpy", "delaunay", "knn_k5"), function(i) {

  grid::grid.grabExpr(pheatmap::pheatmap(cell_proximity_matrix[[i]],
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     angle_col = 45,
                     fontsize = 6,
                     legend=FALSE))
})

p_legends <- lapply(c("squidpy", "delaunay", "knn_k5"), function(i) {
  
  pheatmap::pheatmap(cell_proximity_matrix[[i]],
                                         cluster_rows = FALSE,
                                         cluster_cols = FALSE,
                                         angle_col = 45,
                                         fontsize = 6,
                                         main = i, legend=TRUE)
})

legends <- gridExtra::grid.arrange(p_legends[[1]]$gtable$grobs[[5]],
                        p_legends[[2]]$gtable$grobs[[5]],
                        p_legends[[3]]$gtable$grobs[[5]],
                        ncol=3)
ggsave("plots/vizgen/vizgen_giotto_squidpy_nhood_enrichment_legend.pdf",
       legends,
       width = 6.5, height = 2,
       dpi = 300)

pheatmaps <- cowplot::plot_grid(plotlist = p_heatmaps, ncol=3)

pspatconnect <- cowplot::plot_grid(p_spat_connect_plots[[1]],
                                   p_spat_connect_plots[[2]],
                                   p_spat_connect_plots[[3]],
                                   ncol=3)


both_plots <- cowplot::plot_grid(pspatconnect,
                       pheatmaps,
                       ncol=1
)                   
ggsave("plots/vizgen/vizgen_giotto_squidpy_nhood_enrichment_comparison.pdf",
       both_plots,
       width = 6.5, height = 4,
       dpi = 300)

packageVersion("Giotto")
