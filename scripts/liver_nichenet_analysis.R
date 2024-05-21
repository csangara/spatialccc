# Liver_nichenet
library(Seurat)
library(tidyverse)
library(nichenetr)
library(Giotto)

liver_seurat_obj <- readRDS("LiverData/LiverData_Raw.rds")
liver_seurat_obj <- readRDS("LiverData/LiverData_RawNorm_NormalLiver_Q1.rds")
liver_giotto_obj <- readRDS("LiverData/LiverData_RawNorm_NormalLiver_Q1_giotto.rds")
cell_proximities <- readRDS("cell_proximities_NormalLiver_Q1.rds")

liver_seurat_obj <- liver_seurat_obj[,liver_seurat_obj$Run_Tissue_name == "NormalLiver"]

ggplot(liver_seurat_obj@meta.data %>% filter(Run_Tissue_name != "NormalLiver"),
            aes(x = x_slide_mm, y = y_slide_mm, color = cellType)) +
  geom_point(size=0.05) + coord_fixed() +
  theme_classic()
ggsave("LiverData/CancerousLiver_celltypes.png",
       width = 3000, height = 3000, units = "px")

ggplot(liver_seurat_obj@meta.data %>% filter(Run_Tissue_name == "NormalLiver"),
       aes(x = x_slide_mm, y = y_slide_mm)) +
  geom_point(size=0.05) + facet_wrap(vars(cellType)) +
  theme_classic()

ggsave("NormalLiver_cellType_facetplot.png", units = "mm", width = 297*2, height = 210*2)


p1 <- ggplot(liver_seurat_obj@meta.data %>% mutate(nCount_RNA = as.numeric(nCount_RNA)),
             aes(x=nCount_RNA, color=Run_Tissue_name)) + geom_density() + theme_classic() +
     ggtitle("Count distribution")
p2 <- ggplot(liver_seurat_obj@meta.data %>% mutate(nFeature_RNA = as.numeric(nFeature_RNA)),
             aes(x=nFeature_RNA, color=Run_Tissue_name)) + geom_density() + theme_classic() +
    ggtitle("Feature distribution")
p1 + p2 + patchwork::plot_layout(guides="collect") & theme(legend.position = "bottom",
                                                           axis.title = element_blank()) 
ggsave("count_and_feature_dist.png", width=3000, height=1500, units="px")
# First things first - see overlap of 1000 genes and NN database
ligand_target_matrix <- readRDS("ligand_target_matrix_nsga2r_final.rds")
lr_network <- readRDS("lr_network_human_21122021.rds")

gene_names <- rownames(liver_seurat_obj)
ligand_target_matrix_subset <- ligand_target_matrix[, colnames(ligand_target_matrix) %in% gene_names ]

ligands <- rownames(liver_seurat_obj)[rownames(liver_seurat_obj) %in% (lr_network %>% pull(from))]
receptors <- rownames(liver_seurat_obj)[rownames(liver_seurat_obj) %in% (lr_network %>% pull(to))]

length(ligands) # 319
length(receptors) # 268
length(rownames(liver_seurat_obj)[rownames(liver_seurat_obj) %in% rownames(ligand_target_matrix)]) # all genes in ligand_target_matrix
intersect(ligands, receptors) %>% length # 57
set.seed(20)
random_ligands <- lt_df$ligand %>% unique %>% sample(20)


# lt_df <- data.frame(ligand_target_matrix) %>% rownames_to_column("target") %>% pivot_longer(!target) %>% rename(ligand=name) %>%
#   mutate(scenario = case_when(ligand %in% gene_names & target %in% gene_names ~ "both",
#                               ligand %in% gene_names & !(target %in% gene_names) ~ "only_ligand",
#                               !(ligand %in% gene_names) & target %in% gene_names ~ "only_target",
#                               !(ligand %in% gene_names) & !(target %in% gene_names) ~ "neither"))

lt_df <- data.frame(ligand_target_matrix_subset) %>% rownames_to_column("target") %>% pivot_longer(!target) %>% rename(ligand=name) %>%
  mutate(target_in_network = target %in% gene_names)

# ggplot(lt_df, aes(y=scenario, x=value)) +
#   geom_violin() + theme_classic()

ggplot(lt_df %>% filter(ligand %in% random_ligands), aes(x=ligand, y=value, color=target_in_network)) +
  geom_violin() + theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("random_ligand_plot.png", width=2500, height=1500, units="px")

ggplot(lt_df, aes(x=target_in_network, y=value)) +
  geom_violin() + theme_classic() +
  labs(x="Target gene is in CosMX dataset", title = "Ligand-target matrix filtered to 319 ligands", y = "Regulatory potential")
ggsave("all_ligand_Target_plot.png", width=1500, height=1000, units="px")

grep("BMP", rownames(liver_seurat_obj), value = TRUE)
grep("VSIG4", rownames(liver_seurat_obj), value = TRUE)
grep("GLUL", rownames(liver_seurat_obj), value = TRUE)
grep("EPCAM", rownames(liver_seurat_obj), value = TRUE)

grep("CD", rownames(liver_seurat_obj), value = TRUE)


liver_seurat_obj$gene <- as.vector(GetAssayData(liver_seurat_obj["GLUL",]))

ggplot(liver_seurat_obj@meta.data %>% arrange(gene),
       #aes(x= approximateumap_2, y = approximateumap_1, color = gene)) +
  aes(x = x_slide_mm, y = y_slide_mm, color = gene)) +
  scale_color_gradient(low = "gray90", high = "purple") +
  geom_point(size=0.2) +
  theme_classic()

ggplot(liver_seurat_obj@meta.data ,
       aes(x = x_slide_mm, y = y_slide_mm, color = cellType)) +
  geom_point(size=0.2) +
  theme_classic()

##### INFLAMMATORY VS NON-INFLAMMATORY MACS DE ######
liver_seurat_obj
Idents(liver_seurat_obj)<-liver_seurat_obj$cellType

liver_giotto_obj <- createGiottoObject(GetAssayData(liver_seurat_obj, slot="counts"),
                                       norm_expr = GetAssayData(liver_seurat_obj),
                                       spatial_locs = liver_seurat_obj@meta.data %>% select(x_slide_mm, y_slide_mm),
                                       cell_metadata = liver_seurat_obj@meta.data %>% select(!c(x_slide_mm, y_slide_mm)))
liver_giotto_obj <- createSpatialNetwork(gobject = liver_giotto_obj, minimum_k = 2)
liver_giotto_obj <- createSpatialKNNnetwork(gobject = liver_giotto_obj, k=5)

spatPlot(gobject = liver_giotto_obj,
         show_network = T,
         point_shape = "no_border",
         network_color = 'black',
         spatial_network_name = 'Delaunay_network',
         point_size = 3,
         cell_color = "cellType",
         coord_fix_ratio = 1)

plotStatDelaunayNetwork(gobject = liver_giotto_obj)
cell_proximities <- cellProximityEnrichment(gobject = liver_giotto_obj,
                                           cluster_column = 'cellType',
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

cell_proximities$enrichm_res %>% data.frame %>% filter(grepl("Non.inflammatory", cell_1))

cell_proximities$enrichm_res %>% data.frame %>% filter(grepl("Non.inflammatory", cell_2) | grepl("Non.inflammatory", cell_1)) %>%
  filter(type_int == "hetero", enrichm > 0) %>% arrange(-enrichm) %>% select(cell_1, cell_2) %>% unlist %>% .[!grepl("Non.inflammatory", .)] %>% unname
cell_proximities$enrichm_res %>% data.frame %>% filter(grepl("Inflammatory", cell_2) | grepl("Inflammatory", cell_1)) %>%
  filter(type_int == "hetero", enrichm > 0) %>% arrange(-enrichm) %>% select(cell_1, cell_2) %>% unlist %>% .[!grepl("Inflammatory", .)] %>% unname


# NicheNet
receiver <- "Inflammatory.macrophages"
expressed_genes_receiver <- get_expressed_genes(receiver, liver_seurat_obj, pct = 0.10)
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]


geneset_oi <- FindMarkers(liver_seurat_obj, ident.1 = "Inflammatory.macrophages", ident.2 = "Non.inflammatory.macrophages", slot = "counts")
geneset_oi <- geneset_oi %>% data.frame %>% filter(p_val_adj <= 0.05) %>% rownames %>% .[. %in% rownames(ligand_target_matrix)]

## sender
sender_celltypes <- c("Antibody.secreting.B.cells", "Cholangiocytes", "Periportal.LSECs", "Stellate.cells")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, liver_seurat_obj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
ligands <- lr_network %>% pull(from) %>% unique()
receptors <- lr_network %>% pull(to) %>% unique()

expressed_ligands <- intersect(ligands,expressed_genes_sender)
expressed_receptors <- intersect(receptors,expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

ligand_activities <- predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr) %>% mutate(rank = rank(desc(aupr)))
ligand_activities

# possible_sender_celltypes <- liver_seurat_obj$cellType %>% unique %>% .[!grepl("macrophages", .)]

# all_ligand_activities <- lapply(possible_sender_celltypes, function(sender_celltype) {
#   
#   expressed_genes_sender <- get_expressed_genes(sender_celltype, liver_seurat_obj, 0.10)
#   
#   ligands <- lr_network %>% pull(from) %>% unique()
#   receptors <- lr_network %>% pull(to) %>% unique()
#   
#   expressed_ligands <- intersect(ligands,expressed_genes_sender)
#   expressed_receptors <- intersect(receptors,expressed_genes_receiver)
#   
#   potential_ligands <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
#   
#   ligand_activities <- predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes,
#                                                  ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
#   
#   ligand_activities <- ligand_activities %>% arrange(-aupr) %>% mutate(rank = rank(desc(aupr)))
#   ligand_activities
#   
# })

ggplot(liver_seurat_obj %>% AddMetaData(GetAssayData(liver_seurat_obj)["CSF1",], col.name="ligand") %>% .@meta.data %>%
         arrange(ligand),
       aes(x = x_slide_mm, y = y_slide_mm, color = ligand)) +
  geom_point() + facet_wrap(vars(cellType)) +
  scale_color_gradient(low="gray90", high="red") +
  theme_classic()




best_upstream_ligands = ligand_activities %>% top_n(20, aupr) %>% arrange(-aupr) %>% pull(test_ligand) %>% unique()
DotPlot(liver_seurat_obj, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network

spatGenePlot2D(liver_giotto_obj, genes = c("CSF1", "CSF1R"))

#### COEXPRESSION PLOT ####
liver_giotto_obj <- createSpatialKNNnetwork(gobject = liver_giotto_obj, k=5)
knn_network <- liver_giotto_obj@spatial_network$knn_network$networkDT %>%
                  mutate(from_cellType = liver_seurat_obj@meta.data[.$from, "cellType"],
                         to_cellType = liver_seurat_obj@meta.data[.$to, "cellType"])

knn_network_subset <- knn_network %>% filter(from_cellType == "Non.inflammatory.macrophages" | to_cellType == "Non.inflammatory.macrophages")
  # Create a new column to represent the cell type of the neighboring cells
  # mutate(neighbour_cellType = ifelse(from_cellType == "Non.inflammatory.macrophages", to_cellType, from_cellType))

df_from <- knn_network %>% filter(from_cellType == "Non.inflammatory.macrophages") %>%
  rename(cellType = from_cellType, neighbour_cellType = to_cellType, cell_id = from, neighbor_id = to)
df_to <-  knn_network %>% filter(to_cellType == "Non.inflammatory.macrophages") %>%
  rename(cellType = to_cellType, neighbour_cellType = from_cellType, cell_id = to, neighbor_id = from)

count_df <- bind_rows(df_from, df_to) %>% group_by(cell_id, neighbour_cellType) %>%  summarise(num_neighbours = n())

count_df %>%
  pivot_wider(
    id_cols = cell_id,
    names_from = neighbour_cellType,
    values_from = num_neighbours,
    values_fill = 0
  )

count_df

ggplot(count_df, aes(y=neighbour_cellType, x=num_neighbours)) +
  geom_violin() + theme_classic()

ggplot(liver_seurat_obj %>% AddMetaData(GetAssayData(liver_seurat_obj)["CSF1",], col.name="ligand") %>% .@meta.data %>%
         arrange(ligand),
       aes(x = x_slide_mm, y = y_slide_mm, color = ligand)) +
  geom_point() + facet_wrap(vars(cellType)) +
  scale_color_gradient(low="gray90", high="red") +
  theme_classic()

coexpr_df <- bind_rows(df_from, df_to) %>% mutate("ligand_expr" = t(data.frame(liver_seurat_obj["CSF1",]@assays$RNA@counts))[.$neighbor_id,],
                                                  "target_expr" = t(data.frame(liver_seurat_obj["CSF1R",]@assays$RNA@counts))[.$cell_id,],
                                                  "coexpression" = ligand_expr*target_expr)


# coexpr_df <- bind_rows(df_from, df_to) %>% mutate("ligand_expr" = t(data.frame(liver_seurat_obj["CSF1", .$neighbor_id]@assays$RNA@data)),
#                                      "target_expr" = t(data.frame(liver_seurat_obj["CSF1R", .$cell_id]@assays$RNA@data)),
#                                      "coexpression" = ligand_expr*target_expr)

ggplot(coexpr_df, aes(x=distance, y=coexpression)) + 
  geom_point() +
  facet_wrap(~neighbour_cellType, scales = "free_x") + theme_classic()

ggplot(coexpr_df, aes(x=ligand_expr, y=target_expr)) + 
  geom_point() +
  facet_wrap(~neighbour_cellType, scales = "free") + theme_classic()

ggplot(data.frame(count = liver_seurat_obj$nCount_RNA), aes(x=count)) + geom_density() + theme_classic()
