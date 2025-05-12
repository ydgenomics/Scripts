# Title: metaNeighbor.R
# Date: 2025-05-12
# 基于RNA/SCT做分析

library(MetaNeighbor)
library(SummarizedExperiment)
library(Seurat)
library(SingleCellExperiment)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input_file"),
    type = "character", default = "/data/work/integration/input/Peanut-unsoupx.cg.rds",
    help = "Path to input file"
  ),
  make_option(c("-o", "--output_name"),
    type = "character", default = "peanut",
    help = "Output file prefix name"
  ),
  make_option(c("-b", "--batch_key"),
    type = "character", default = "biosample",
    help = "Batch key for integration"
  ),
  make_option(c("-c", "--cluster_key"),
    type = "character", default = "leiden_res_0.50",
    help = "Cluster key for integration"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))
input_file <- opt$input_file
out_put_name <- opt$output_name
batch_key <- opt$batch_key
cluster_key <- opt$cluster_key


sdata <- readRDS(input_file)
sdata
colnames(sdata@meta.data)
sdata <- as.SingleCellExperiment(sdata)
head(colData(sdata))

var_genes = variableGenes(dat = sdata, exp_labels = sdata@colData[[batch_key]])

celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = sdata,
                             study_id = sdata@colData[[batch_key]],
                             cell_type = sdata@colData[[cluster_key]],
                             fast_version = TRUE)


cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

#in work directory output a pdf
pdf(paste0(out_put_name,"_metaNeighbor.pdf"))
gplots::heatmap.2(celltype_NV,
                  margins=c(8,8),
                  keysize=1,
                  key.xlab="AUROC",
                  key.title="NULL",
                  trace = "none",
                  density.info = "none",
                  col = cols,
                  breaks = breaks,
                  offsetRow=0.1,
                  offsetCol=0.1,
                  cexRow = 0.7,
                  cexCol = 0.7)
dev.off()

top_hits = topHits(cell_NV = celltype_NV,
                   dat = sdata,
                   study_id = sdata$species,
                   cell_type = sdata$celltype,
                   threshold = 0.9)

write.csv(file=paste0(out_put_name,"_metaNeighbor_tophits.csv"),top_hits,quote=FALSE,row.names=FALSE)