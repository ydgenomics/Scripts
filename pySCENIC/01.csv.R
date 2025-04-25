# Title: 01.csv.R
# Author: ydgenomics
# Date: 2025-04-22
# Description: This script is the preparation of runing pySCENIC.
# Pre-requiry:
# params: rds_path

library(Seurat)
library(stringr)

rds_path <- 

seu <- readRDS(rds_path)
colnames(seu@meta.data)
head(t(as.matrix(seu@assays$RNA@counts))[1:3,1:3])
head(colnames(t(as.matrix(seu@assays$RNA@counts))))

# save csv
write.csv(t(as.matrix(seu@assays$RNA@counts)),file = "scenic.data.csv")
