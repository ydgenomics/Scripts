# Introduction
# DecontPro assess and decontaminate single-cell protein expression data, such as 
# those generated from CITE-seq or Total-seq. The count matrix is decomposed into 
# three matrices, the native, the ambient and the background that represent the 
# contribution from the true protein expression on cells, the ambient material and
# other non-specific background contamination.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("decontX")
library(decontX)