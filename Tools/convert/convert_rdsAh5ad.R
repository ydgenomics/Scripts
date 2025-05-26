# Title: convert_rdsAh5ad.R
# Date: 20250526
# Coder: ydgenomics
# Description: Using sceasy(R) and schardto convert single-cell data between rds and h5ad
# Image: /software/conda/Anaconda/bin/R
# Reference: Â© EMBL-European Bioinformatics Institute, 2023 Yuyao Song <ysong@ebi.ac.uk>


library(sceasy)
library(reticulate)
library(Seurat)
library(schard)
packageVersion("Seurat")
library(optparse)

option_list <- list(
  make_option(c("-i", "--input_file"),type = "character", default = NULL, help = "Path to input file for convrting"),
  make_option(c("-a", "--assay"),type = "character", default = "RNA",help = "Assay name for the output file"),
  make_option(c("-m", "--main_layer"),type = "character", default = "counts",help = "Main layer name for the output file")
)
opt <- parse_args(OptionParser(option_list = option_list))
input_path <- opt$input_file
assay <- opt$assay
main_layer <- opt$main_layer

convert_rdsAh5ad <- function(input_path, assay = "RNA", main_layer = "counts") {
    ext <- tools::file_ext(input_path)
    print(paste0("input file extension is : ", ext))
    if (ext == "rds") {
        file_name <- basename(input_path)
        output_path <- sub("\\.rds$", ".rh.h5ad", file_name)
        message(paste0("from seurat to anndata, input: ", input_path))
        use_python("/opt/conda/bin/python")
        loompy <- reticulate::import("loompy")
        temp0 <- readRDS(input_path)
        print(temp0)
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
        temp0[[assay]] <- as(temp0[[assay]], "Assay")
        sceasy::convertFormat(
            temp0,
            from = "seurat",
            to = "anndata",
            assay = assay,
            main_layer = main_layer,
            outFile = output_path
        )
        cat("Converted rds to h5ad. Output: ", output_path, "\n")
        return(output_path)
    } else if (ext == "h5ad") {
        file_name <- basename(input_path)
        output_path <- sub("\\.h5ad$", ".hr.rds", file_name)
        message(paste0("from anndata to seurat, input: ", input_path))
        dt <- schard::h5ad2seurat(input_path)
        saveRDS(dt, file = output_path)
        cat("Converted h5ad to rds. Output: ", output_path, "\n")
        return(output_path)
    } else {
        stop("Error: The file extension is neither .rds nor .h5ad.")
    }
}

output_file <- tryCatch({
  convert_rdsAh5ad(input_path, assay, main_layer)
}, error = function(e) {
  cat("An error occurred: ", e$message, ", but don't worried. \n")
  NULL
}, warning = function(w) {
  cat("A warning occurred: ", w$message, "\n")
  NULL
}, finally = {
  cat("Attempt to convert RDS to H5AD completed.\n")
  NULL
})