# Title: convert_rdsAh5ad2.R
# Date: 20250529
# Coder: ydgenomics
# Description: Using sceasy(R) and schard to convert single-cell data between Seurat and AnnData formats.
# Image: /software/conda/Anaconda/bin/R
# Reference: Â© EMBL-European Bioinformatics Institute, 2023 Yuyao Song <ysong@ebi.ac.uk>


library(sceasy)
library(reticulate)
library(Seurat)
library(schard)
packageVersion("Seurat")
library(optparse)

option_list <- list(
  make_option(
    c("-i", "--input_file"),type = "character", 
    default = "/data/work/0.peanut/annotation/three_layers/H1314_dataget_Anno_rename_threelayers.h5ad", 
    help = "Path to input file for convrting"),
  make_option(
    c("-l", "--layers"), type = "character",
    default = "RNA",
    help = "Layers to be converted, default is RNA")
)
opt <- parse_args(OptionParser(option_list = option_list))
input_path <- opt$input_file
layers <- opt$layers
ext <- tools::file_ext(input_path)
print(paste0("input file extension is : ", ext))
layers <- unlist(strsplit(layers, ",")); print(paste0("layers to be converted: ", paste(layers, collapse = ",")))


if (ext == "rds") {
    message(paste0("from seurat to anndata, input: ", input_path))
    # Using reticulate to call Python
    use_python("/opt/conda/bin/python")
    loompy <- reticulate::import("loompy")
    temp0 <- readRDS(input_path)
    # Check '_index' in meta.data and meta.features
    if ("_index" %in% colnames(temp0@meta.data)) {
        print("Change _index into new_index for meta.data")
        colnames(temp0@meta.data)[colnames(temp0@meta.data) == "_index"] <- "new_index"
    }
    if ("meta.features" %in% slotNames(temp0[["RNA"]])) {
        features_meta <- temp0[["RNA"]]@meta.features
        if ("_index" %in% colnames(features_meta)) {
            colnames(features_meta)[colnames(features_meta) == "_index"] <- "new_index"
            print("Change _index into new_index for meta.features")
        }
        temp0[["RNA"]]@meta.features <- features_meta
    } else {
        message("meta.features does not exist in the RNA Assay of temp0.")
    }
    colnames(temp0@meta.data)
    assays <- names(temp0@assays)
    if (layers[1] == "all"){
        print("Layers is set to 'all', using all assays in the Seurat object.")
    } else {
        assays <- intersect(assays, layers)
        print(paste0("Layers is set to: ", paste(assays, collapse = ",")))
        if (length(assays) == 0) {
            stop("Error: No matching layers found in the Seurat object.")
        }
    }
    h5ad_paths <- c()
    file_name <- basename(input_path)
    output_path <- sub("\\.rds$", ".rh.h5ad", file_name)
    file_name_without_ext <- sub("\\.rds$", "", file_name)
    getwdpath <- getwd()
    for (i in assays) {
        h5ad_path <- sprintf("%s_%s_rh.h5ad", file_name_without_ext, i)
        h5ad_path <- paste0(getwdpath, "/", h5ad_path)
        sceasy::convertFormat(
            temp0,
            from = "seurat",
            to = "anndata",
            assay = i,
            main_layer = "counts",
            outFile = h5ad_path
        )
        cat("Converted rds to h5ad. Output: ", h5ad_path, "\n")
        h5ad_paths <- c(h5ad_paths, h5ad_path)
    }
    writeLines(paste(h5ad_paths, collapse = ","), "saved_paths.txt")
    writeLines(paste(assays, collapse = ","), "saved_layers.txt")
} else if (ext == "h5ad") {
    file_name <- basename(input_path)
    output_path <- sub("\\.h5ad$", ".hr.rds", file_name)
    message(paste0("from anndata to seurat, input: ", input_path))
    source("/data/work/multi_anno/AT_root/convert_rdsAh5ad.R")
    # 调用 Python 函数
    saved_layers <- unlist(strsplit(readLines("saved_layers.txt"), ","))
    if ("counts" %in% saved_layers) {saved_layers[saved_layers == "counts"] <- "RNA"}; print(saved_layers)
    saved_paths <- unlist(strsplit(readLines("saved_paths.txt"), ",")); print(saved_paths)
    if (layers[1] == "all"){
        print("Layers is set to 'all', using all assays in the Seurat object.")
        layers <- saved_layers
    } else {
        print(paste0("Layers is set to: ", paste(layers, collapse = ",")))
    }
    common_layers <- intersect(saved_layers, layers)
    filtered_saved_layers <- saved_layers[saved_layers %in% common_layers]
    filtered_saved_paths <- saved_paths[saved_layers %in% common_layers]
    print("Filtered saved_layers:"); print(filtered_saved_layers)
    print("Filtered saved_paths:"); print(filtered_saved_paths)
    rds_paths <- c()
    for (path in filtered_saved_paths) {
        cat("Processing file:", path, "\n")
        rds_path <- convert_rdsAh5ad(path)
        rds_paths <- c(rds_paths, rds_path)
    }
    print(rds_paths)
    seu <- readRDS(rds_paths[1])
    if (length(filtered_saved_paths) < 2) {
        print(paste0("Only one layer is present: ", filtered_saved_paths[1]))
    } else {
        print(paste0("Multiple layers are present: ", paste(filtered_saved_paths, collapse = ", ")))
        for (i in 2:length(filtered_saved_paths)) {
            seu2 <- readRDS(rds_paths[i])
            saved_layer <- saved_layers[i]
            rna_data <- GetAssayData(seu2, assay = "RNA", layer = "counts")
            other_assay <- CreateAssayObject(counts = rna_data, meta.data = seu2@meta.data, name = saved_layer)
            seu[[saved_layer]] <- other_assay
            cat(sprintf("Layer: %s, Added to Seurat object\n", saved_layer))
        }   
    }
    print(seu)
    DefaultAssay(seu) <- "RNA"
    saveRDS(seu, output_path)
    # Delete the temporary files
    files_to_delete <- c(saved_paths, rds_paths)
    for (file_path in files_to_delete) {
      if (file.exists(file_path)) {
        file.remove(file_path)
        cat("Deleted file:", file_path, "\n")
      } else {
        cat("File does not exist, could not delete:", file_path, "\n")
      }
    }
    cat("Converted h5ad to rds. Output: ", output_path, "\n")
} else {
    stop("Error: The file extension is neither .rds nor .h5ad.")
}