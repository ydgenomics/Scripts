library(Seurat)
input_query_rds <- "/data/users/yangdong/yangdong_f6fd22e6e3a247ceaae97934225564ba/online/result/scPlant/test/NipLSD1_obj_after_choir.rds"
input_marker_csv <- "/data/work/multi_anno/rice_leaf_marker.csv"
input_ref_rdata <- "/data/users/yangdong/yangdong_f6fd22e6e3a247ceaae97934225564ba/online/result/scPlant/Reference_CRA004082_leaf_filtered_cg_pp.Rdata"
seu <- readRDS(input_query_rds)
source("/data/work/multi_anno/run_sctype.R")
seu <- run_sctype(seu, cluster_key="CHOIR_clusters_0.05", input_marker_csv, tissue="leaf", umap_name="CHOIR_P0_reduction_UMAP")

source("/data/work/multi_anno/run_singler.R")
seu <- run_singler(seu, input_ref_rdata)
