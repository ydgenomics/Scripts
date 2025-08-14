# Date: 250814
# Image: metaNeighbor-R--04 /opt/conda/bin/R
# Reference: https://mp.weixin.qq.com/s/tVxalBWsxLn58RJkpb-PaQ
# 基于RNA@counts做分析

library(MetaNeighbor)
library(SummarizedExperiment)
library(Seurat)
library(SingleCellExperiment)
library(ambient)
library(grid)
library(ComplexHeatmap)
library(ggcor)
library(circlize)
library(ggplot2)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input_file"),
    type = "character", default = "/data/work/SCPipelines/bulk_RNA_scRNA_singleR/bulk_anno_seu.rds",
    help = "Path to input file"
  ),
  make_option(c("-o", "--output_name"),
    type = "character", default = "cotton",
    help = "Output file prefix name"
  ),
  make_option(c("-b", "--batch_key"),
    type = "character", default = "time",
    help = "Batch key for integration"
  ),
  make_option(c("-c", "--cluster_key"),
    type = "character", default = "leiden_res_0.50",
    help = "Cluster key for integration"
  ),
  make_option(c("-t", "--threshold_value"),
    type = "numeric", default = 0.95,
    help = "Threshold value for top hits"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))
input_file <- opt$input_file
output_name <- opt$output_name
batch_key <- opt$batch_key
cluster_key <- opt$cluster_key
threshold_value <- opt$threshold_value


sdata <- readRDS(input_file)
sdata
colnames(sdata@meta.data)
sdata <- as.SingleCellExperiment(sdata, assay = "RNA", slot = "counts")
# print(assay(sdata, "counts")[1:5, 1:5])
head(colData(sdata))

var_genes = variableGenes(dat = sdata, exp_labels = sdata@colData[[batch_key]])

celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = sdata,
                             study_id = sdata@colData[[batch_key]],
                             cell_type = sdata@colData[[cluster_key]],
                             fast_version = TRUE)

# Using ggcor do circle heatmap
p1 <- quickcor(celltype_NV, circular = TRUE, cluster = TRUE, grid.colour = 'white',
         open = 90, # 缺口大小
         # 内圈外圈比例
         outer = 0.1, inner = 0.2) +
  # 单元格边框线颜色
  geom_colour(colour = 'black') +
  # 自定义填充颜色
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') +
  # 更改图例名称
  guides(fill = guide_colorbar(title = 'AUROC')) +
  anno_col_tree() +
  anno_row_tree() +
  # 基因名
  set_p_yaxis() +
  # 样本名
  set_p_xaxis()

p2 <- quickcor(celltype_NV, circular = TRUE, cluster = TRUE,
         open = 90, # 缺口大小
         # 内圈外圈比例
         outer = 0.1, inner = 0.3) +
  # 单元格边框线颜色
  geom_colour(colour = 'black') +
  # 自定义填充颜色
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') +
  # 更改图例名称
  guides(fill = guide_colorbar(title = 'AUROC')) +
  # 添加聚类树
  anno_col_tree(height = 0.05, bcols = c('#0A81AB','#FB9300')) +
  anno_row_tree(pos = 'left', bcols = rainbow(5)) +
  # 基因名
  set_p_yaxis() +
  # 样本名
  set_p_xaxis()

p3 <- quickcor(celltype_NV, circular = TRUE, cluster = TRUE, grid.colour = 'white',
         open = 90, # 缺口大小
         # 内圈外圈比例
         outer = 0.2, inner = 0.3) +
  # 单元格边框线颜色
  geom_colour(colour = 'white') +
  # 自定义填充颜色
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') +
  # 更改图例名称
  guides(fill = guide_colorbar(title = 'AUROC')) +
  # 列注释
  anno_hc_bar(k = 2, fill = rand_color(2), pos = 'top', height = 0.3) +
  anno_hc_bar(k = 3, fill = rand_color(3), pos = 'top', height = 0.5) +
  anno_hc_bar(k = 5, fill = rand_color(5), pos = 'top', height = 0.2) +
  # 添加聚类树
  anno_col_tree(bcols = rand_color(5), height = 0.15) +
  anno_hc_bar(k = 15, fill = rand_color(15), pos = 'left', width = 0.5) +
  anno_hc_bar(k = 10, fill = rand_color(10), pos = 'left', width = 0.5) +
  anno_row_tree(bcols = rand_color(8)) +

  # 样本名
  set_p_xaxis(bcols = rand_color(5)) +
  # 基因名
  set_p_yaxis()


cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

# in work directory output a pdf
pdf(paste0(output_name,"_metaNeighbor.pdf"))
# Using heatmap do heatmap
gplots::heatmap.2(celltype_NV,
                  col = cols,
                  breaks = breaks,
                  key.xlab = "AUROC",
                  margins = c(8, 8),
                  trace = "none",
                  density.info = "none",
                  offsetRow=0.1,
                  offsetCol=0.1,
                  cexRow = 0.7,
                  cexCol = 0.7)

# ggcor
print(p1)
print(p2)
print(p3)

dev.off()
write.csv(file=paste0(output_name,"_metaNeighbor.csv"),celltype_NV,quote=FALSE,row.names=TRUE)
top_hits = topHits(cell_NV = celltype_NV,
                   dat = sdata,
                   study_id = sdata@colData[[batch_key]],
                   cell_type = sdata@colData[[cluster_key]],
                   threshold = threshold_value)

write.csv(file=paste0(output_name,"_metaNeighbor_tophits.csv"),top_hits,quote=FALSE,row.names=FALSE)