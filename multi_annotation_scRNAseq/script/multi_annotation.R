# parameters
marker_csv_path="/data/work/multi_anno/rice_leaf_marker.csv"
tissue="leaf"
umap_name="CHOIR_P0_reduction_UMAP"
cluster_key="CHOIR_clusters_0.05"
input_query_rds <- "/data/users/yangdong/yangdong_f6fd22e6e3a247ceaae97934225564ba/online/result/scPlant/test/NipLSD1_obj_after_choir.rds"
input_ref_rdata <- "/data/users/yangdong/yangdong_f6fd22e6e3a247ceaae97934225564ba/online/result/scPlant/Reference_CRA004082_leaf_filtered_cg_pp.Rdata"

# load rds file
seu <- readRDS(input_query_rds); DefaultAssay(seu) <- 'RNA'

# preprocess
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)

# run sctype
source("/data/work/multi_anno/run_sctype.R")
seu <- run_sctype(seu, cluster_key=cluster_key, input_marker_csv=marker_csv_path, tissue=tissue, umap_name=umap_name)

# run SingleR
source("/data/work/multi_anno/run_singler.R")
load(input_ref_rdata) # load reference data
seu <- run_singler(seu, ref_sce)

# save result
colnames(seu@meta.data)
output_query_rds <- paste0(sub("\\.rds$", "", basename(input_query_rds)), "_sctype.rds")
#output_umap <- paste0(sub("\\.rds$", "", basename(input_query_rds)), "_sctype_umap.pdf")

output_umap <- "output_sctype_umap.pdf"
pdf(output_umap)
DimPlot(seu, reduction = umap_name, label = TRUE, pt.size = 0.5, group.by = cluster_key)
DimPlot(seu, reduction = umap_name, label = TRUE, pt.size = 0.5, group.by = "sctype")
DimPlot(seu, reduction = umap_name, label = TRUE, pt.size = 0.5, group.by = "singler")
dev.off()
saveRDS(seu, file = "/data/work/multi_anno/seu_multi_annotation.rds")
