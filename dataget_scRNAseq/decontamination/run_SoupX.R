# Date: 20250606
# Image: SoupX-R--03 /opt/conda/bin/R
# Reference: https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html

library(optparse)
library(DropletUtils)
library(SoupX)
library(Seurat)

option_list <- list(
    make_option(c("-r", "--raw_path"), type = "character", default = "/data/input/Files/cuiyingsi/my-result2/v3/V3RNA25012000009/YS2-V3RNA25012000009/output/raw_matrix", help = "String: Path to raw matrix", metavar = "character"),
    make_option(c("-f", "--filter_path"), type = "character", default = "/data/input/Files/cuiyingsi/my-result2/v3/V3RNA25012000009/YS2-V3RNA25012000009/output/filter_matrix", help = "String: Path to filtered matrix", metavar = "character"),
    make_option(c("-s", "--sample_name"), type = "character", default = "V3RNA25012000009", help = "String: Sample name", metavar = "character"),
    make_option(c("-m", "--minCG"), type = "numeric", default = 100, help = "Minimum number of genes", metavar = "numeric"),
    make_option(c("-t", "--tfidfMin"), type = "numeric", default = 0.01, help = "Minimum tf-idf value", metavar = "numeric"),
    make_option(c("-x", "--highestrho"), type = "numeric", default = 0.2, help = "Highest acceptable rho value", metavar = "numeric")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# raw_path <- opt$raw_path
# filter_path <- opt$filter_path
# sample_name <- opt$sample_name
# minCG <- opt$minCG
# tfidfMin <- opt$tfidfMin
# highestrho <- opt$highestrho
raw_path <- "/data/input/Files/cuiyingsi/my-result2/v3/V3RNA25012000010/YS2-V3RNA25012000010/output/raw_matrix"
filter_path <- "/data/input/Files/cuiyingsi/my-result2/v3/V3RNA25012000010/YS2-V3RNA25012000010/output/filter_matrix"
sample_name <- "V3RNA25012000010"
minCG <- 100
tfidfMin <- 0.01
highestrho <- 0.2

run_soupx <- function(raw_path, filter_path, sample_name, minCG, tfidfMin, highestrho) {
    options(future.globals.maxSize = 100000 * 1024^3)
    toc <- Read10X(filter_path, gene.column=1)
    # Because CreateSeuratObject() will replace '_' as '-', in order to keep raw genes' name
    gene_names <- rownames(toc); head(gene_names)
    all <- CreateSeuratObject(toc)
    rownames(all) <- gene_names; head(rownames(all))
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

    tod <- Read10X(raw_path, gene.column=1)
    gene_names <- rownames(tod); head(gene_names)
    raw <- CreateSeuratObject(tod)
    rownames(raw) <- gene_names; head(rownames(raw))
    tod <- GetAssayData(object = raw, layer = "counts", assay = "RNA")
    tod <- tod[rownames(all),]

    sc <- SoupChannel(tod, toc)
    #sc <- estimateSoup(sc)
    sc <- setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
    pdfrho <- paste0(sample_name, "_rho.pdf")
    pdf(pdfrho, width=9)
    sc <- autoEstCont(sc, tfidfMin = tfidfMin, forceAccept = TRUE) # If you want mannually set rho, use setContaminationFraction(), but it is not recommended
    dev.off()

    rho_value <- unique(sc$metaData$rho)
    rho_adjust <- ifelse(rho_value < highestrho, 'good', 'bad')

    out <- adjustCounts(sc)
    DropletUtils::write10xCounts(sample_name, out, version="3")
    file_conn <- file(paste0(sample_name,'_soupx_rho.txt'), open = "w")
    cat("highestrho:", highestrho, "\n", file = file_conn)
    cat("sample:", sample_name, "\n", file = file_conn)
    cat("rho:", rho_value, "\n", file = file_conn) 
    cat("how:", rho_adjust, "\n", file = file_conn)
    close(file_conn)
}

result <- run_soupx(raw_path=raw_path, filter_path=filter_path, sample_name=sample_name, minCG=minCG, tfidfMin=tfidfMin, highestrho=highestrho)