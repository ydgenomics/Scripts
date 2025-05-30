input_marker_csv="/data/work/multi_anno/rice_leaf_marker.csv"
tissue="leaf"
umap_name="CHOIR_P0_reduction_UMAP"
input_query_rds <- "/data/users/yangdong/yangdong_f6fd22e6e3a247ceaae97934225564ba/online/result/scPlant/test/NipLSD1_obj_after_choir.rds"

seu <- readRDS(input_query_rds); DefaultAssay(seu) <- 'RNA'
source("/data/work/multi_anno/run_sctype.R")
seu <- run_sctype(seu, cluster_key="CHOIR_clusters_0.05", input_marker_csv=input_marker_csv, tissue=tissue, umap_name=umap_name)
input_ref_rdata <- "/data/users/yangdong/yangdong_f6fd22e6e3a247ceaae97934225564ba/online/result/scPlant/Reference_CRA004082_leaf_filtered_cg_pp.Rdata"
load(input_ref_rdata)
source("/data/work/multi_anno/run_singler.R")
seu <- run_singler(seu, ref_sce)
colnames(seu@meta.data)
