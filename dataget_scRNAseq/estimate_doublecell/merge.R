### Date: 250725
### Image: sceasy-schard--
### Coder: ydgenomics

library(Seurat)
library(dplyr)
# library(presto)
library(optparse)

option_list <- list(
    make_option(c("-p", "--prefix"), type = "character", default = "test",
                help = "Prefix for output files [default %default]"),
    make_option(c("-r", "--rds_files"), type = "character",
                default = "/data/input/W202507240000066_1/Dataget/call-merge/execution/peanut_merge/H1314.hr.rds,/data/input/W202507240000066_1/Dataget/call-merge/execution/peanut_merge/H2014.hr.rds",
                help = "Comma-separated list of RDS files [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
prefix <- opt$prefix
rds_files <- opt$rds_files

### Read RDS files
rds_files <- unlist(strsplit(rds_files, ",")); print(rds_files)
if (length(rds_files) > 1) {
    seu <- readRDS(rds_files[[1]])
    for (i in 2:length(rds_files)) {
        temp_data <- readRDS(rds_files[[i]])
        print(head(colnames(temp_data)))
        seu <- merge(seu, temp_data)
    }
}
print("------------------------- Merged Seurat Object ----------------------")
print(seu)
print("---- meta.data columns ----")
print(colnames(seu@meta.data))
print("---- meta.data head ----")
print(head(seu@meta.data))

print("--------------- Normalization and Clustering -----------------")
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, nfeatures = 3000)
# seu <- FindVariableFeatures(seu, selection.method = "mean.var.plot", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, features = VariableFeatures(object = seu), verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5, cluster.name = "merge_res0.5")
print("------------- meta.data columns after clustering -------------")
print(colnames(seu@meta.data))
seu <- RunUMAP(seu, dims = 1:20, verbose = FALSE)
#
variable_features <- VariableFeatures(seu)
print("---- Variable Features head ----")
print(head(variable_features))

pdf(paste0(prefix, "_merge.pdf"), width = 10, height = 8)
DimPlot(seu, reduction = "umap", group.by = "biosample", shuffle = TRUE, label = TRUE)
DimPlot(seu, reduction = "umap", group.by = "sample", shuffle = TRUE, label = TRUE)
# 'leiden_res_0.20', 'leiden_res_0.50', 'leiden_res_0.80', 'leiden_res_1.00', 'leiden_res_1.30', 'leiden_res_1.60', 'leiden_res_2.00'
DimPlot(seu, reduction = "umap", group.by = "leiden_res_0.50", shuffle = TRUE, label = TRUE)
DimPlot(seu, reduction = "umap", group.by = "leiden_res_0.80", shuffle = TRUE, label = TRUE)
DimPlot(seu, reduction = "umap", group.by = "leiden_res_1.00", shuffle = TRUE, label = TRUE)
DimPlot(seu, reduction = "umap", group.by = "leiden_res_1.30", shuffle = TRUE, label = TRUE)
DimPlot(seu, reduction = "umap", group.by = "leiden_res_1.60", shuffle = TRUE, label = TRUE)
DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", shuffle = TRUE, label = TRUE)
dev.off()

saveRDS(seu, paste0(prefix,"_merge.rds"))