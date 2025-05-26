# Title: convert_rdsAh5ad2.R
# Date: 20250526
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
  make_option(c("-i", "--input_file"),type = "character", 
  default = "/data/work/convert/0525test/multi_layers/H1314_dataget_Anno_rename_threelayers.hr.rds", 
  help = "Path to input file for convrting")
)
opt <- parse_args(OptionParser(option_list = option_list))
input_path <- opt$input_file
ext <- tools::file_ext(input_path)
print(paste0("input file extension is : ", ext))

# Using reticulate to call Python functions
use_python("/opt/conda/bin/python")
source_python("/script/deal_layers_ydgenomics.py")

if (ext == "rds") {
    message(paste0("from seurat to anndata, input: ", input_path))
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
    h5ad_paths <- c()
    file_name <- basename(input_path)
    output_path <- sub("\\.rds$", ".rh.h5ad", file_name)
    file_name_without_ext <- sub("\\.rds$", "", file_name)
    for (i in assays) {
        h5ad_path <- sprintf("%s_%s_rh.h5ad", file_name_without_ext, i)
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
    py$get_multi_layers_h5ad(h5ad_paths, assays, output_path)
    # 遍历每个文件路径并删除文件
    for (file_path in h5ad_paths) {
      if (file.exists(file_path)) {
        file.remove(file_path)
        cat("Deleted file:", file_path, "\n")
      } else {
        cat("File does not exist, could not delete:", file_path, "\n")
      }
    }
} else if (ext == "h5ad") {
    file_name <- basename(input_path)
    output_path <- sub("\\.h5ad$", ".hr.rds", file_name)
    message(paste0("from anndata to seurat, input: ", input_path))
    source("/script/convert_rdsAh5ad.R")
    # 调用 Python 函数
    layers_list <- py$get_single_layer_h5ad(input_path)
    saved_layers <- layers_list[[1]]
    saved_paths <- layers_list[[2]]
    rds_paths <- c()
    for (path in saved_paths) {
        cat("Processing file:", path, "\n")
        rds_path <- convert_rdsAh5ad(path)
        rds_paths <- c(rds_paths, rds_path)
    }
    print(rds_paths)
    seu <- readRDS(rds_paths[1])
    for (i in 2:length(saved_layers)) {
        seu2 <- readRDS(rds_paths[i])
        saved_layer <- saved_layers[i]
        rna_data <- GetAssayData(seu2, assay = "RNA", layer = "counts")
        other_assay <- CreateAssayObject(counts = rna_data, meta.data = seu2@meta.data, name = saved_layer)
        seu[[saved_layer]] <- other_assay
        cat(sprintf("Layer: %s, Added to Seurat object\n", saved_layer))
    }
    print(seu)
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