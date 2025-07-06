# convert_format.R 250706
# /software/conda/Anaconda/bin/R
# Â© EMBL-European Bioinformatics Institute, 2023
# Yuyao Song <ysong@ebi.ac.uk>

library(Seurat)
packageVersion("Seurat")
library(schard)
library(sceasy)
library(reticulate)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input_file"),
    type = "character", default = "/data/work/Arabidopsis.hr.rds",
    help = "Path to input file for convrting"
  ),
  make_option(c("-o", "--output_file"),
    type = "character", default = "/data/work/Arabidopsis.hr.rh.h5ad",
    help = "Output file after conversion"
  ),
  make_option(c("-a", "--assay"),
    type = "character", default = "RNA",
    help = "Assay name for the output file"
  ),
  make_option(c("-m", "--main_layer"),
    type = "character", default = "counts",
    help = "Main layer name for the output file"
  ),
  make_option(c("-s", "--stype"),
    type = "character", default = "rds2h5ad",
    help = "One format convert to another format"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
input_file <- opt$input_file
output_file <- opt$output_file
assay <- opt$assay
main_layer <- opt$main_layer
stype <- opt$stype


if(stype == 'h5ad2rds'){
    message(paste0("from anndata to seurat, input: ", input_file))
    dt=schard::h5ad2seurat(input_file)
    saveRDS(dt,file=output_file)
} else if (stype == 'rds2h5ad'){
    message(paste0("from seurat to anndata, input: ", input_file))
    use_python("/software/conda/Anaconda/bin/python")
    loompy <- reticulate::import('loompy')
    temp0 <- readRDS(input_file)
    DefaultAssay(temp0) <- 'RNA'
    # Check '_index' in meta.data and meta.feature, because RLIGER way will product '_index' cause error of python dataframe
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
    message("Meta.data info of inputed seurat object")
    print(colnames(temp0@meta.data))
    temp0[["RNA"]] <- as(temp0[["RNA"]], "Assay")
    sceasy::convertFormat(temp0, from="seurat", to="anndata", assay = assay, main_layer=main_layer, outFile = output_file)
}