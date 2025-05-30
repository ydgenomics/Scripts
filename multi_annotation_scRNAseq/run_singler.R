# Title: singler.R
# Date: 20250528
# Coder: ydgenomics
# Description: Using SingleR to annotate single-cell RNA-seq data based on a custom reference dataset.
# Input: reference .rds has RNA, query .rds files, and a metadata key for clustering in the reference dataset
# Output: `_refsingleR.Rdata` and `_singleR.rds files`
# Image: Seurat-R--04 /software/miniconda/envs/Seurat/bin/R
# Reference: [使用singleR基于自建数据库来自动化注释单细胞转录组亚群](https://mp.weixin.qq.com/s/GpOxe4WLIrBOjbdH5gfyOQ)

library(Seurat)
library(SingleCellExperiment)
library(scater)
library(SingleR)
library(optparse)

# option_list <- list(
#     make_option(
#         c("-r", "--input_ref_rds"), type = "character", default = "data/immune_ref.rds", help = "Path to the reference dataset"),
#     make_option(
#         c("-q", "--input_query_rds"), type = "character", default = "data/immune_query.rds", help = "Path to the query dataset"),
#     make_option(
#         c("-k", "--ref_cluster_key"), type = "character", default = "Celltype",help = "Metadata key for clustering in the reference dataset")
# )
# opt <- parse_args(OptionParser(option_list = option_list))
# input_ref_rds <- opt$input_ref_rds
# input_query_rds <- opt$input_query_rds
# ref_cluster_key <- opt$ref_cluster_key

# Step 1: Load the reference dataset and create a singleR reference Rdata object
create_ref_singler <- function(ref_seu, ref_cluster_key) {
    #ref_seu <- readRDS(input_ref_rds)
    print(ref_seu); print(colnames(ref_seu))
    Idents(ref_seu) <- ref_seu@meta.data[[ref_cluster_key]]
    av <- AggregateExpression(
        ref_seu, group.by = ref_cluster_key, assays = "RNA"
    )
    ref_mat <- av[[1]]
    ref_sce <- SingleCellExperiment(
        assays = list(counts = ref_mat)
    )
    ref_sce <- scater::logNormCounts(ref_sce)
    colData(ref_sce)$Type <- colnames(ref_mat)
    output_ref_rdata <- paste0(sub("\\.rds$", "", basename(input_ref_rds)), "_refsingleR.Rdata")
    save(ref_sce, file = output_ref_rdata)
    return(ref_sce)
}


# Step 2: Load the query dataset and run singleR for annotation
run_singler <- function(query_seu, ref_sce) {
    #query_seu <- readRDS(input_query_rds); DefaultAssay(query_seu) <- "RNA"
    DefaultAssay(query_seu) <- "RNA"
    query_seu <- NormalizeData(query_seu, normalization.method = "LogNormalize", scale.factor = 10000)
    query_data <- GetAssayData(query_seu, slot = "data")
    common_genes <- intersect(rownames(query_data), rownames(ref_sce))
    num_common_genes <- length(common_genes)
    print(
        paste0(
            "The common gene number of Reference and Query data: ",
            num_common_genes
        )
    )
    pred <- SingleR(
        test = query_data, ref = ref_sce, labels = ref_sce$Type
    )
    query_seu$singleR <- pred$labels
    #output_query_rds <- paste0(sub("\\.rds$", "", basename(input_query_rds)), "_singleR.rds")
    #saveRDS(query_seu, file = output_query_rds)
    return(query_seu)
}

# run_singler(input_query_rds, ref_sce)