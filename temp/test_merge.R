# Date: 250628
# 1.探究merge的用法，以及对比python的concat函数的异同
# 2.探究先NormalizeData再merge和先merge再NormalizeData的区别
# 3.实践如何基于2个表创建seurat对象

setwd("/data/work/script/merge0627")
library(Seurat)
library(dplyr)

name <- "merge0627"
file_paths <- c("/data/work/script/rename0626/D1/D1_0626.rds",
               "/data/work/script/rename0626/D2/D2_0626.rds",
               "/data/work/script/rename0626/D3/D3_0626.rds",
               "/data/work/script/rename0626/D5/D5_0626.rds",
               "/data/work/script/rename0626/D10/D10_0626.rds",
               "/data/work/script/rename0626/D30/D30_0626.rds")

if (length(file_paths) > 1) {
    merged_data <- readRDS(file_paths[[1]])
    merged_data$time <- names[[1]]
    print(head(colnames(merged_data)))
    for (i in 2:length(file_paths)) {
        temp_data <- readRDS(file_paths[[i]])
        print(paste0("Processing: ", names[[i]]))
        temp_data$time <- names[[i]]
        print(head(colnames(temp_data)))
        merged_data <- merge(merged_data, temp_data)
    }
}

head(merged_data$RNA@counts)
head(merged_data$RNA@data)
colnames(merged_data@meta.data)
#  [1] "orig.ident"          "nCount_RNA"          "nFeature_RNA"       
#  [4] "percent.mt"          "DF.classifications"  "RNA_snn_res.0.5"    
#  [7] "seurat_clusters"     "cell_type"           "cell_type_marker"   
# [10] "sample"              "new_clusters"        "CHOIR_clusters_0.05"
# [13] "sctype"              "rename0626"          "time"               
# [16] "L16_clusters"        "L19_clusters"       
unique(merged_data$time)
unique(merged_data$orig.ident)
#################### 直接使用merge前的data信息 ##################
seu <- merged_data
seu
seu <- FindVariableFeatures(seu, nfeatures = 3000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, features = VariableFeatures(object = seu), verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5, cluster.name = "merge0627_res0.5")
seu <- RunUMAP(seu, dims = 1:20, verbose = FALSE)
variable_features <- VariableFeatures(seu)
print(head(variable_features))
pdf(paste0(name, "_before_normalize_merged_umap.pdf"), width = 16, height = 8)
DimPlot(seu, reduction = "umap", group.by = "time", shuffle = TRUE, label = TRUE)
DimPlot(seu, reduction = "umap", group.by = "sample", shuffle = TRUE, label = TRUE)
DimPlot(seu, reduction = "umap", group.by = "rename0626", shuffle = TRUE, label = TRUE)
DimPlot(seu, reduction = "umap", group.by = "merge0627_res0.5", shuffle = TRUE, label = TRUE)
dev.off()
saveRDS(seu, paste0(name,"_before_normalize_merged.rds"))

#################### 一起NormalizeData ##################
merged_data <- NormalizeData(merged_data)
merged_data <- FindVariableFeatures(merged_data, nfeatures = 3000)
merged_data <- ScaleData(merged_data)
merged_data <- RunPCA(merged_data, features = VariableFeatures(object = merged_data), verbose = FALSE)
merged_data <- FindNeighbors(merged_data, dims = 1:30)
merged_data <- FindClusters(merged_data, resolution = 0.5, cluster.name = "merge0627_res0.5")
merged_data <- RunUMAP(merged_data, dims = 1:20, verbose = FALSE)
#
variable_features <- VariableFeatures(merged_data)
print(head(variable_features))
pdf(paste0(name, "_after_normalize_merged_umap.pdf"), width = 16, height = 8)
DimPlot(merged_data, reduction = "umap", group.by = "time", shuffle = TRUE, label = TRUE)
DimPlot(merged_data, reduction = "umap", group.by = "sample", shuffle = TRUE, label = TRUE)
DimPlot(merged_data, reduction = "umap", group.by = "rename0626", shuffle = TRUE, label = TRUE)
DimPlot(merged_data, reduction = "umap", group.by = "merge0627_res0.5", shuffle = TRUE, label = TRUE)
dev.off()

saveRDS(merged_data, paste0(name,"_after_normalize_merged.rds"))