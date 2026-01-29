library(tidyverse)
library(patchwork)
library(Seurat)
library(sf)

leiden_cluster_palette <- c(RColorBrewer::brewer.pal(n = 12, name = "Set3"))

color_palette <- c("Hepatocytes" = "#B4B5B5FF",
                   "Endothelial cells" = "#dca754",
                   "Cholangiocytes" = "#C61B84FF",
                   "Stromal cells" = "#79151e",
                   "Kuppfer cells" = "#5DA6DBFF",
                   "Other immunecells" = "#893A86FF",
                   "B cells" = "#9C7EBAFF",
                   "Unknown" = "#191919")


# Read vizgen data
vizgen_data <- readRDS("data/vizgen/vizgen_data.rds")
vizgen_data@meta.data %>% head
GetAssayData(vizgen_data)[1:5,1:5]
GetAssayData(vizgen_data, layer = "counts")[1:5,1:5]

DimPlot(vizgen_data, group.by = "leiden")
DimPlot(vizgen_data, group.by = "annotation_own_score_genes")

leiden_other_res <- read.csv("rds/vizgen_leiden_clusters.csv")
vizgen_data$leiden_0.3 <- as.factor(leiden_other_res$leiden_0.3)
vizgen_data$leiden_0.4 <- as.factor(leiden_other_res$leiden_0.4)
vizgen_data$leiden_0.5 <- as.factor(leiden_other_res$leiden_0.5)
vizgen_data$celltype <- vizgen_data$annotation_own_score_genes

# Read niche detection results
nichecompass <- read.csv("rds/vizgen_nichecompass_latent_clusters_and_umap.csv")
covet <- read.csv("rds/vizgen_covet_clusters_and_umap.csv")

all(Cells(vizgen_data) == nichecompass$cells)
all(Cells(vizgen_data) == covet$cells)

vizgen_data$nichecompass_cluster <- as.factor(nichecompass$latent_leiden_0.4)
vizgen_data$covet_cluster <- as.factor(covet$leiden_covet)

DimPlot(vizgen_data, group.by = "nichecompass_cluster")
DimPlot(vizgen_data, group.by = "covet_cluster")

# Plot ROI
roi <- c(10000,15000,8000,12000)
cols_oi <- c("celltype", "leiden_0.3", "nichecompass_cluster", "covet_cluster")

segmented_cells <- read_sf("rds/vizgen_segmentation_mask_boundaries_roi9-16-7-13.shp")

p_spatial <- lapply(cols_oi, function(tool) {
  
  if (tool == "celltype"){
    palette <- color_palette
  } else {
    palette <- leiden_cluster_palette
  }
  
  selected_cluster <- vizgen_data@meta.data[tool] %>% 
    rename(selected_cluster = all_of(tool)) %>%
    rownames_to_column("index")
  
  segmented_cells_join <- segmented_cells %>%
    left_join(selected_cluster, by="index")
  
  ggplot() +
    geom_sf(aes(fill=selected_cluster), color = "gray50", linewidth=0.1, data = segmented_cells_join) +
    theme_void(base_size=7) +
    scale_fill_manual(values=palette, name = "Leiden cluster", drop=FALSE) +
    labs(title=tool) +
    theme(#legend.position = "none",
      plot.title = element_text(hjust=0.5, size=8),
      legend.key.size = unit(0.3, "cm"))
  
}) %>% setNames(cols_oi)

# Create UMAPs


umap_meta <- bind_cols(Embeddings(vizgen_data, reduction = "umap"),
          vizgen_data@meta.data[,cols_oi])
umap_meta_long <- umap_meta %>%
  rownames_to_column("cell") %>%
  pivot_longer(cols = all_of(cols_oi), names_to = "tool", values_to = "cluster") 

umap_median <- umap_meta_long %>% 
  group_by(tool, cluster) %>% 
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

p_umap_cluster <- lapply(cols_oi, function(tool_i) {
  
  if (tool_i == "celltype"){
    palette <- color_palette
  } else {
    palette <- leiden_cluster_palette
  }
  
  ggplot(umap_meta_long %>% filter(tool == tool_i)) +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = cluster), 
             size=0.2, stroke=0, shape=16) +
  # geom_label(data = umap_median %>% filter(tool == tool_i),
  #            aes(x = UMAP_1, y = UMAP_2, label = cluster, fill = cluster),
  #            color = "black", size=2) +
  guides(color = guide_legend(override.aes = list(size=1.5))) +
  scale_color_manual(values = palette) +
  #scale_fill_manual(values = palette) +
  ggtitle(paste0("Leiden clusters")) +
  theme_classic(base_size=7) +
  theme(plot.title = element_text(hjust = 0.5, size=7, face="bold"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")
  
}) %>% setNames(cols_oi)

wrap_plots(p_umap_cluster)
# Create boxplots

props_df <- vizgen_data@meta.data %>% select(all_of(cols_oi)) %>% 
  rownames_to_column("cell") %>% 
  rename(celltype = annotation_own_score_genes) %>%
  pivot_longer(-c(cell,celltype), names_to = "tool", values_to = "cluster") %>% 
  group_by(tool, celltype, cluster) %>% 
  summarise(n=n()) %>% 
  group_by(tool, cluster) %>% 
  mutate(freq = n/sum(n))

p_barplots <- lapply(unique(props_df$tool), function(tool_i){
  ggplot(props_df %>% filter(tool == tool_i),
         aes(y=forcats::fct_rev(cluster), x=n, fill=celltype)) +
    geom_bar(position = "stack", stat = "identity", width=0.5) +
    theme_minimal(base_size=7) +
    scale_fill_manual(values=color_palette, name="Cell type") +
    ggtitle(tool_i) +
    labs(x="Frequency", y="Leiden clusters") +
    theme(legend.key.size = unit(0.3, "cm"),
          legend.position = "bottom",
          axis.text.x = element_text(angle=0),
          panel.grid.major.y = element_blank())

}) %>% setNames(unique(props_df$tool))

p_spatial[[1]]
wrap_plots(p_spatial)
wrap_plots(p_barplots)

cols_oi

design <- "JK#
           ABC
           DEF
           GHI"

p_spatial[["leiden_0.3"]] + p_spatial[["nichecompass_cluster"]] + p_spatial[["covet_cluster"]]
p_umap_cluster[["leiden_0.3"]] + p_umap_cluster[["nichecompass_cluster"]] + p_umap_cluster[["covet_cluster"]] &
  theme(legend.position = "right")
p_umap_cluster[["celltype"]]
wrap_plots(p_barplots)

p_final <- p_spatial[["leiden_0.3"]] + p_spatial[["nichecompass_cluster"]] + p_spatial[["covet_cluster"]] +
  p_umap_cluster[["leiden_0.3"]] + p_umap_cluster[["nichecompass_cluster"]] + p_umap_cluster[["covet_cluster"]] +
  p_barplots[["leiden_0.3"]] + p_barplots[["nichecompass_cluster"]] + p_barplots[["covet_cluster"]] +
  p_umap_cluster[["celltype"]] + p_spatial[["celltype"]] +
  plot_layout(design = design, guides = "collect") &
  theme(legend.position = "none", plot.title = element_blank(),
        axis.title = element_blank())
p_final

ggsave("plots/vizgen/vizgen_niche_detection_comparison_2.pdf",
       plot = p_final,
       width = 6.5, height = 7,
       dpi = 300)


# Get legend
barplot_legend <- cowplot::get_legend(
  p_barplots[[1]] +
    guides(fill = guide_legend(title = "Cell type", nrow=1)) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size=5, margin = margin(l=2, r=4)),
          legend.title = element_text(size=6))
)
cowplot::plot_grid(barplot_legend)



spatial_legend_celltype <- 
  cowplot::get_legend(
  p_spatial[["celltype"]] +
    guides(fill = guide_legend(title = "Celltype",
                               override.aes = list(linewidth = 0, color="black"))) +
    theme(legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size=5),
          legend.title = element_text(size=6))
    )

spatial_legend_leiden <- 
  cowplot::get_legend(
    p_spatial[["leiden_0.3"]] +
      guides(fill = guide_legend(title = "Cluster", nrow=1,
                                 override.aes = list(linewidth = 0, color="black"))) +
      theme(legend.position = "bottom",
            legend.direction = "horizontal",
            legend.text = element_text(size=5, margin = margin(l=2, r=4)),
            legend.key.size = unit(0.3, "cm"),
            legend.title = element_text(size=6))
  )
spatial_legend_grid <- cowplot::plot_grid(spatial_legend_celltype,
                                          spatial_legend_leiden, nrow=1)

spatial_legend_grid
# Plot legends with cowplot
legends <- cowplot::plot_grid(
  barplot_legend,
  spatial_legend_grid,
  nrow=2
)

legends

ggsave("plots/vizgen/vizgen_niche_detection_comparison_legend.pdf",
       plot = legends,
       width = 6.5, height = 4,
       dpi = 300)
