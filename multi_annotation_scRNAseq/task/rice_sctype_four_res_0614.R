# parameters
marker_csv_path="/data/work/script/sctype0614/rice_leaf_marker0614.csv"
tissue="leaf"
umap_name="CHOIR_P0_reduction_UMAP"
cluster_key="CHOIR_clusters_0.05"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript rice_sctype_four_res_0614.R <input_query_rds>")
}
input_query_rds <- args[1]

# run sctype
source("/data/work/script/sctype0614/run_sctype.R")

# load rds file
seu <- readRDS(input_query_rds); DefaultAssay(seu) <- 'RNA'
seu <- FindClusters(seu, resolution = 0.5, cluster.name = "resolution0.5"); unique(seu$resolution0.5)
seu <- FindClusters(seu, resolution = 0.8, cluster.name = "resolution0.8"); unique(seu$resolution0.8)
seu <- FindClusters(seu, resolution = 1.0, cluster.name = "resolution1.0"); unique(seu$resolution1.0)
# # preprocess
# seu <- NormalizeData(seu)
# seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
# seu <- ScaleData(seu)
# obj <- FindClusters(obj, resolution = resolution_set, cluster.name = "celltype")


cluster_key="CHOIR_clusters_0.05"
seu <- run_sctype(seu, cluster_key=cluster_key, input_marker_csv=marker_csv_path, tissue=tissue, umap_name=umap_name)
anno_key <- paste0("sctype_",cluster_key)
seu@meta.data[[anno_key]] <- seu$sctype
anno_key2 <- paste0(anno_key, "_", cluster_key)
seu@meta.data[[anno_key2]] <- paste0(seu@meta.data[[anno_key]],"_",seu@meta.data[[cluster_key]])
colnames(seu@meta.data)
file.rename("output_sctype_umap.pdf", paste0(sub("\\.rds$", "", basename(input_query_rds)), "_", cluster_key,"_sctype_umap.pdf"))

cluster_key <- "resolution0.5"
seu <- run_sctype(seu, cluster_key=cluster_key, input_marker_csv=marker_csv_path, tissue=tissue, umap_name=umap_name)
anno_key <- paste0("sctype_",cluster_key)
seu@meta.data[[anno_key]] <- seu$sctype
anno_key2 <- paste0(anno_key, "_", cluster_key)
seu@meta.data[[anno_key2]] <- paste0(seu@meta.data[[anno_key]],"_",seu@meta.data[[cluster_key]])
colnames(seu@meta.data)
file.rename("output_sctype_umap.pdf", paste0(sub("\\.rds$", "", basename(input_query_rds)), "_", cluster_key,"_sctype_umap.pdf"))

cluster_key <- "resolution0.8"
seu <- run_sctype(seu, cluster_key=cluster_key, input_marker_csv=marker_csv_path, tissue=tissue, umap_name=umap_name)
anno_key <- paste0("sctype_",cluster_key)
seu@meta.data[[anno_key]] <- seu$sctype
anno_key2 <- paste0(anno_key, "_", cluster_key)
seu@meta.data[[anno_key2]] <- paste0(seu@meta.data[[anno_key]],"_",seu@meta.data[[cluster_key]])
colnames(seu@meta.data)
file.rename("output_sctype_umap.pdf", paste0(sub("\\.rds$", "", basename(input_query_rds)), "_", cluster_key,"_sctype_umap.pdf"))

cluster_key <- "resolution1.0"
seu <- run_sctype(seu, cluster_key=cluster_key, input_marker_csv=marker_csv_path, tissue=tissue, umap_name=umap_name)
anno_key <- paste0("sctype_",cluster_key)
seu@meta.data[[anno_key]] <- seu$sctype
anno_key2 <- paste0(anno_key, "_", cluster_key)
seu@meta.data[[anno_key2]] <- paste0(seu@meta.data[[anno_key]],"_",seu@meta.data[[cluster_key]])
colnames(seu@meta.data)
file.rename("output_sctype_umap.pdf", paste0(sub("\\.rds$", "", basename(input_query_rds)), "_", cluster_key,"_sctype_umap.pdf"))

output_umap <- paste0(sub("\\.rds$", "", basename(input_query_rds)), "_sctype_umap.pdf")
pdf(output_umap, width=14, height=8)
DimPlot(seu, reduction = umap_name, label = TRUE, pt.size = 0.5, group.by = "sctype_CHOIR_clusters_0.05")
DimPlot(seu, reduction = umap_name, label = TRUE, pt.size = 0.5, group.by = "sctype_CHOIR_clusters_0.05_CHOIR_clusters_0.05")
DimPlot(seu, reduction = umap_name, label = TRUE, pt.size = 0.5, group.by = "sctype_resolution0.5")
DimPlot(seu, reduction = umap_name, label = TRUE, pt.size = 0.5, group.by = "sctype_resolution0.5_resolution0.5")
DimPlot(seu, reduction = umap_name, label = TRUE, pt.size = 0.5, group.by = "sctype_resolution0.8")
DimPlot(seu, reduction = umap_name, label = TRUE, pt.size = 0.5, group.by = "sctype_resolution0.8_resolution0.8")
DimPlot(seu, reduction = umap_name, label = TRUE, pt.size = 0.5, group.by = "sctype_resolution1.0")
DimPlot(seu, reduction = umap_name, label = TRUE, pt.size = 0.5, group.by = "sctype_resolution1.0_resolution1.0")
dev.off()

colnames(seu@meta.data)
output_query_rds <- paste0(sub("\\.rds$", "", basename(input_query_rds)), "_sctype.rds")
saveRDS(seu, file = output_query_rds)