# Title: 01.csv.R
# Author: ydgenomics
# Date: 20250520
# Description: This script is the preparation of runing pySCENIC.
# Pre-requiry:
# params: rds_path

library(Seurat)
library(stringr)
library(optparse)

option_list <- list(
    make_option(c("-r", "--rds"), type = "character", default = "input.rds",help = "Path to the input RDS file", metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
rds <- opt$rds

seu <- readRDS(rds)
colnames(seu@meta.data)
head(t(as.matrix(seu@assays$RNA@counts))[1:3,1:3])
head(colnames(t(as.matrix(seu@assays$RNA@counts))))

# save csv
write.csv(t(as.matrix(seu@assays$RNA@counts)),file = "scenic.data.csv")
