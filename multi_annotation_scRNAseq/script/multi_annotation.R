# parameters
marker_csv_path="/data/work/script/sctype0614/rice_leaf_marker0614.csv"
tissue="leaf"
umap_name="CHOIR_P0_reduction_UMAP"
cluster_key="CHOIR_clusters_0.05"
input_query_rds <- "/data/input/Files/huangpeilin/NipLSD-result-new1/1/NipLSD1_anno_merged_data_obj_after_choir.rds"
input_ref_rdata <- "/data/users/yangdong/yangdong_f6fd22e6e3a247ceaae97934225564ba/online/result/scPlant/Reference_CRA004082_leaf_filtered_cg_pp.Rdata"

library(Seurat)
# load rds file
seu <- readRDS(input_query_rds); DefaultAssay(seu) <- 'RNA'
print(seu)

# # preprocess
# seu <- NormalizeData(seu)
# seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
# seu <- ScaleData(seu)

# run sctype
source("/data/work/test0615/run_sctype.R")
seu <- run_sctype(seu, cluster_key=cluster_key, input_marker_csv=marker_csv_path, tissue=tissue, umap_name=umap_name)
anno_key1 <- paste0("sctype_", cluster_key)
seu@meta.data[[anno_key1]] <- paste0(seu$sctype,"_",seu@meta.data[[cluster_key]])

# run SingleR
source("/data/work/test0615/run_singler.R")
load(input_ref_rdata) # load reference data
seu <- run_singler(seu, ref_sce)
anno_key2 <- paste0("singler_", cluster_key)
seu@meta.data[[anno_key2]] <- paste0(seu$sctype,"_",seu@meta.data[[cluster_key]])

# save result
colnames(seu@meta.data)

output_umap <- paste0(sub("\\.rds$", "", basename(input_query_rds)), "_sctype_umap.pdf")
pdf(output_umap, width=14, height=8)
DimPlot(seu, reduction = umap_name, label = TRUE, pt.size = 0.5, group.by = cluster_key)
DimPlot(seu, reduction = umap_name, label = TRUE, pt.size = 0.5, group.by = "sctype")
DimPlot(seu, reduction = umap_name, label = TRUE, pt.size = 0.5, group.by = anno_key1)
DimPlot(seu, reduction = umap_name, label = TRUE, pt.size = 0.5, group.by = "singler")
DimPlot(seu, reduction = umap_name, label = TRUE, pt.size = 0.5, group.by = anno_key2)
dev.off()
output_query_rds <- paste0(sub("\\.rds$", "", basename(input_query_rds)), "_sctype.rds")
saveRDS(seu, file = output_query_rds)