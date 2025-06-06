# Title: soupx.R 
# Date: 2025-04-14
# /opt/conda/bin/R
library(optparse)
library(DropletUtils)
library(SoupX)
library(Seurat)

option_list <- list(
    make_option(c("-r", "--raw"), type = "character", default = "raw_matrix.txt", help = "Path to raw matrix file", metavar = "character"),
    make_option(c("-f", "--filter"), type = "character", default = "filter_matrix.txt", help = "Path to filtered matrix file", metavar = "character"),
    make_option(c("-s", "--samples"), type = "character", default = "samples.txt", help = "Path to samples file", metavar = "character"),
    make_option(c("-m", "--minCG"), type = "numeric", default = 200, help = "Minimum number of genes", metavar = "numeric"),
    make_option(c("-t", "--tfidfMin"), type = "numeric", default = 0.01, help = "Minimum tf-idf value", metavar = "numeric"),
    make_option(c("-x", "--highestrho"), type = "numeric", default = 0.2, help = "Highest acceptable rho value", metavar = "numeric")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
raw_filepath_txt <- opt$raw
filter_filepath_txt <- opt$filter
samples_txt_path <- opt$samples
minCG <- opt$minCG
tfidfMin <- opt$tfidfMin
highestrho <- opt$highestrho

# Raw matrix
rawmatrix_list <- readLines(raw_filepath_txt, warn = FALSE)
rawmatrix_list <- unlist(strsplit(rawmatrix_list, ',', fixed = TRUE))
print(rawmatrix_list)

# Filter Matrix
filtermatrix_list <- readLines(filter_filepath_txt, warn = FALSE)
filtermatrix_list <- unlist(strsplit(filtermatrix_list, ',', fixed = TRUE))
print(filtermatrix_list)

# Sample names
sample_list <- readLines(samples_txt_path, warn = FALSE)
sample_list <- unlist(strsplit(sample_list, ',', fixed = TRUE))
print(sample_list)

run_soupx <- function(rawMatrix, filterMatrix, outdir, minCG, tfidfMin, highestrho) {
    options(future.globals.maxSize = 100000 * 1024^3)
    tod <- Read10X(rawMatrix, gene.column=1)
    toc <- Read10X(filterMatrix, gene.column=1)
    
    all <- CreateSeuratObject(toc)
    all <- subset(all, subset = nFeature_RNA > minCG)
    all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
    all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
    all.genes <- rownames(all)
    all <- ScaleData(all, features = all.genes)
    all <- RunPCA(all, features = VariableFeatures(all), npcs = 40, verbose = FALSE)
    all <- FindNeighbors(all, dims = 1:30)
    all <- FindClusters(all, resolution = 0.5)
    all <- RunUMAP(all, dims = 1:30)
    
    matx <- all@meta.data
    toc <- GetAssayData(object = all, layer = "counts", assay = "RNA")
    
    raw <- CreateSeuratObject(tod)
    tod <- GetAssayData(object = raw, layer = "counts", assay = "RNA")
    tod <- tod[rownames(all),]
    
    sc <- SoupChannel(tod, toc)
    sc <- setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
    pdfrho <- paste0(outdir, "_rho.pdf")
    pdf(pdfrho, width=9)
    sc <- autoEstCont(sc, tfidfMin = tfidfMin, forceAccept = TRUE)
    dev.off()
    
    rho_value <- unique(sc$metaData$rho)
    rho_adjust <- ifelse(rho_value < highestrho, 'good', 'bad')
    
    out <- adjustCounts(sc)
    DropletUtils::write10xCounts(outdir, out, version="3")
    
    return(list(rho_value = rho_value, rho_adjust = rho_adjust))
}

rho_value_list <- list()
rho_adjust_list <- list()

for (i in 1:length(sample_list)) {
    result <- run_soupx(rawmatrix_list[i], filtermatrix_list[i], sample_list[i], minCG=minCG, tfidfMin=tfidfMin, highestrho=highestrho)
    rho_value_list[[i]] <- result$rho_value
    rho_adjust_list[[i]] <- result$rho_adjust
}

rho_value_list
rho_adjust_list

file_conn <- file('soupx_rho.txt', open = "w")
cat("highestrho:", highestrho, "\n", file = file_conn)
for (i in seq_along(sample_list)) {
    cat("sample:", sample_list[i], "\n", file = file_conn)
    cat("rho:", rho_value_list[[i]][1], "\n", file = file_conn)  # 提取列表中的向量的第一个值
    cat("how:", rho_adjust_list[[i]][1], "\n", file = file_conn)  # 提取列表中的向量的第一个值
}
close(file_conn)