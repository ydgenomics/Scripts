# Date: 20250714
# Image: metaNeighbor-R--02 /opt/conda/bin/R
# Reference: https://mp.weixin.qq.com/s/tVxalBWsxLn58RJkpb-PaQ
# 基于RNA/SCT做分析

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
pdf(paste0(out_put_name,"_metaNeighbor.pdf"))
# Using heatmap do heatmap
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
# ggcor
print(p1)
print(p2)
print(p3)

# # Using circlize do circle-heatmap
# mycol <- colorRamp2(c(0, 0.5, 1), c("#393781", "white", "#f22942"))
# bordercol <- "white"
# # 调整圆环首尾间的距离
# circos.par(gap.after = c(90)) 
# # 绘制环状热图
# circos.heatmap(celltype_NV, col = mycol,
#                # 聚类放在环形内侧
#                dend.side = "inside", 
#                # 基因名放在环形外侧；二者不能在同一侧
#                rownames.side = "outside",
#                rownames.col = "black",
#                # 字体大小
#                rownames.cex = 0.8, 
#                # 字体粗细
#                rownames.font = 2,  # 注意：字体粗细应该是一个整数，例如 1（普通）、2（加粗）
#                bg.border = bordercol,
#                cluster = TRUE,
#                cell.border = bordercol,
#                track.height = 0.5)

# # 定义图例的位置
# x_pos <- 0.75
# y_pos <- 0.65
# # 创建图例
# lg <- Legend(title = "AUROC", col_fun = mycol,
#              direction = "vertical",
#              title_position = "topcenter")

# # 绘制图例
# draw(lg, x = unit(x_pos, "npc"), y = unit(y_pos, "npc"), just = c("right", "center"))

# # 设置间距因子
# spacing_factor <- 0.72
# # 设置 y 轴偏移量，可以根据需要调整
# y_offset <- 2
# # 设置 x 轴偏移量，可以根据需要调整
# x_offset <- -0.3
# fontsize <- 0.6 # 字体大小

# # 添加自定义文本
# circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
#   if (CELL_META$sector.numeric.index == 1) { 
#     cn <- colnames(celltype_NV)
#     n <- length(cn)
#     base_x <- CELL_META$cell.xlim[2] - convert_x(x_offset, "mm") 
#     circos.text(rep(base_x, n), 
#                 y_offset + (n:1) * spacing_factor,
#                 cn, cex = fontsize, adj = c(0, 0.5), facing = "inside")
#   }
# }, bg.border = NA)
# circos.clear()
dev.off()

top_hits = topHits(cell_NV = celltype_NV,
                   dat = sdata,
                   study_id = sdata@colData[[batch_key]],
                   cell_type = sdata@colData[[cluster_key]],
                   threshold = 0.9)

write.csv(file=paste0(out_put_name,"_metaNeighbor_tophits.csv"),top_hits,quote=FALSE,row.names=FALSE)