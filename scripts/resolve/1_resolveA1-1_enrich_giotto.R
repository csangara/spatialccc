library(tidyverse)
library(Giotto)

data_path <- '~/Documents/other_scripts/resolveA1-1/'

metadata <- read.csv(paste0(data_path, "obs.csv"), row.names = 1)
spatial_locs <- read.csv(paste0(data_path, "obsm.csv")) %>% select(spatial1, spatial2) %>% data.frame %>% 
  magrittr::set_rownames(rownames(metadata))
expr <- read.csv(paste0(data_path, "X.csv"), header=FALSE) %>% as.matrix %>% 
  magrittr::set_rownames(rownames(metadata)) %>% 
  magrittr::set_colnames(read.csv(paste0(data_path, "var.csv")) %>% pull(gene))

liver_giotto_obj <- createGiottoObject(expr %>% t,
                                       spatial_locs = spatial_locs,
                                       cell_metadata = metadata)


liver_giotto_obj <- createSpatialNetwork(gobject = liver_giotto_obj, minimum_k = 2)
liver_giotto_obj <- createSpatialKNNnetwork(gobject = liver_giotto_obj, k=5)

spatPlot(gobject = liver_giotto_obj,
         show_network = T,
         point_shape = "no_border",
         network_color = 'black',
         spatial_network_name = 'Delaunay_network',
         point_size = 3,
         cell_color = "annotationSave",
         coord_fix_ratio = 1)

plotStatDelaunayNetwork(gobject = liver_giotto_obj)

cell_proximities2 <- cellProximityEnrichment(gobject = liver_giotto_obj,
                                            cluster_column = 'annotation_own_score_genes',
                                            spatial_network_name = 'Delaunay_network',
                                            adjust_method = 'fdr',
                                            number_of_simulations = 1000)

cell_proximities <- readRDS("cell_proximities_NormalLiver_Q1.rds")
# saveRDS(liver_giotto_obj, "LiverData/LiverData_RawNorm_NormalLiver_Q1_giotto.rds")
# saveRDS(cell_proximities, "cell_proximities_NormalLiver_Q1.rds")

cellProximityBarplot(gobject = liver_giotto_obj,
                     CPscore = cell_proximities,
                     min_orig_ints = 3,
                     min_sim_ints = 3)

cellProximityNetwork(gobject = liver_giotto_obj,
                     CPscore = cell_proximities,
                     remove_self_edges = F,
                     self_loop_strength = 0.3,
                     only_show_enrichment_edges = F,
                     rescale_edge_weights = T,
                     node_size = 8,
                     edge_weight_range_depletion = c(1,2),
                     edge_weight_range_enrichment = c(2,5))

cellProximityHeatmap(gobject = liver_giotto_obj,
                     CPscore = cell_proximities
                     
)

cell_proximities$enrichm_res %>% filter(grepl("Kupffer cells", unified_int))
cell_proximities$enrichm_res %>% filter(grepl("LAM", unified_int))


enrich_res <- cell_proximities$enrichm_res
first_type <- second_type <- unified_int <- NULL
enrich_res[, `:=`(first_type, strsplit(x = as.character(unified_int), 
                                       split = "--")[[1]][1]), by = 1:nrow(enrich_res)]
enrich_res[, `:=`(second_type, strsplit(x = as.character(unified_int), 
                                        split = "--")[[1]][2]), by = 1:nrow(enrich_res)]
enrich_mat <- data.table::dcast.data.table(data = enrich_res, 
                                          formula = first_type ~ second_type, value.var = "enrichm")
matrix_d <- as.matrix(enrich_mat[, -1])
rownames(matrix_d) <- as.vector(enrich_mat[[1]])
t_matrix_d <- Giotto:::t_giotto(matrix_d)
t_matrix_d[upper.tri(t_matrix_d)][is.na(t_matrix_d[upper.tri(t_matrix_d)])] = matrix_d[upper.tri(matrix_d)][is.na(t_matrix_d[upper.tri(t_matrix_d)])]
t_matrix_d[lower.tri(t_matrix_d)][is.na(t_matrix_d[lower.tri(t_matrix_d)])] = matrix_d[lower.tri(matrix_d)][is.na(t_matrix_d[lower.tri(t_matrix_d)])]
t_matrix_d[is.na(t_matrix_d)] = 0
final_matrix = t_matrix_d
final_matrix <-  Giotto:::t_giotto(scale( Giotto:::t_giotto(final_matrix)))
final_matrix <-  Giotto:::t_giotto(final_matrix)
final_matrix[lower.tri(final_matrix)] <-  Giotto:::t_giotto(final_matrix)[lower.tri(final_matrix)]
cordist = stats::as.dist(1 -  Giotto:::cor_giotto(final_matrix))
clus = stats::hclust(cordist)
myorder = clus$order
mylabels = clus$labels
names(mylabels) = 1:length(mylabels)
sample_order = mylabels[myorder]
final_matrix = final_matrix[sample_order, sample_order]


cutree(clus, k=3) %>% sort() %>%
  #split
  split(., .)
  