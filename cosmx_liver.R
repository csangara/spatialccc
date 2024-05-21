library(Seurat)
library(SeuratDisk)

# Scan areas: 76mm2 (normal liver) and 100mm2 (carcinoma)

liver_cosmx_seuratobj <- readRDS("LiverDataReleaseSeurat_newUMAP.RDS")
liver_cosmx_seuratobj <- readRDS("~/spotless-benchmark/unit-test/test_sc_data.rds")


file_name <- stringr::str_split(basename(test[1]), "\\.")[[1]][1]
ext <- stringr::str_split(basename(test[1]), "\\.")[[1]][2]


file_name <- c(Sys.getenv("VSC_DATA_VO_USER"), "/spatialnichenet/LiverDataReleaseSeurat_newUMAP.h5Seurat")

SaveH5Seurat(liver_cosmx_seuratobj, file_name)
