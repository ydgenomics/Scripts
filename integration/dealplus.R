# dealplus.R
# /opt/conda/bin/R
library(Seurat)
library(optparse)
library(ggplot2)

option_list <- list(
  make_option(c("-i", "--input_rds"),
    type = "character", default = NULL,
    help = "Path to input preprocessed rds file"
  ),
  make_option(c("-o", "--out_rds"),
    type = "character", default = NULL,
    help = "integrated rds file"
  ),
  make_option(c("-p", "--out_UMAP"),
    type = "character", default = NULL,
    help = "Output UMAP after integration"
  ),
  make_option(c("-t", "--other1_key"),
    type = "character", default = NULL,
    help = "Output UMAP after integration"
  ),
  make_option(c("-f", "--other2_key"),
    type = "character", default = NULL,
    help = "Output UMAP after integration"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_rds)){
  opt$input_rds <- "/data/work/sct/Cer_test_convert_SCTransform.CCA_integrated.rds"
}
if (is.null(opt$out_rds)){
  opt$out_rds <- "/data/work/dealplus/Cer_test_convert_SCTransform.CCA_integrated.rds"
}
if (is.null(opt$out_UMAP)){
  opt$out_UMAP <- "/data/work/dealplus/Cer_test_convert_SCTransform.CCA_integrated_UMAP.pdf"
}

input_rds <- opt$input_rds
out_rds <- opt$out_rds
out_UMAP <- opt$out_UMAP
other1_key <- opt$other1_key
other2_key <- opt$other2_key

obj <- readRDS(input_rds)
#add a column
obj$biosample.celltype <- paste0(obj$biosample, "_", obj$celltype)
unique(obj$biosample.celltype)

#visual
pdf(out_UMAP, width=12, height=8)
DimPlot(obj, reduction = "umap", split.by = "biosample")
DimPlot(obj, reduction = "umap", group.by = "biosample", shuffle = TRUE, label = TRUE)
DimPlot(obj, reduction = "umap", group.by = "sample", shuffle = TRUE, label = TRUE)
DimPlot(obj, reduction = "umap", group.by = "celltype", shuffle = TRUE, label = TRUE)
DimPlot(obj, reduction = "umap", group.by = "biosample.celltype", shuffle = TRUE, label = TRUE)
# 检查 meta.data 中是否包含指定列
required_features <- c("nCount_RNA", "nFeature_RNA")
if (all(required_features %in% colnames(obj@meta.data))) {
  VlnPlot(obj, features = required_features, group.by = "biosample")
} else {
  message("meta.data lacked nCount_RNA or nFeature_RNA, so passed vioplot")
}

# 检查 meta.data 中是否存在这两个列
if (other1_key %in% colnames(obj@meta.data) && other2_key %in% colnames(obj@meta.data)) {
  # 检查是否满足条件：other1_key 不等于 "biosample" 且 other2_key 不等于 "celltype"
  if (other1_key != "biosample" || other2_key != "celltype") {
    obj$other1_other2 <- paste0(
      obj@meta.data[[other1_key]], "_", obj@meta.data[[other2_key]]
    )
    p1 <- DimPlot(obj, reduction = "umap", split.by = other1_key)
    p2 <- DimPlot(obj, reduction = "umap", group.by = other1_key, shuffle = TRUE, label = TRUE)
    p3 <- DimPlot(obj, reduction = "umap", group.by = other2_key, shuffle = TRUE, label = TRUE)
    p4 <- DimPlot(obj, reduction = "umap", group.by = "other1_other2", shuffle = TRUE, label = TRUE)
    print(p1); print(p2); print(p3); print(p4)
  } else {
    message("other1_key 或 other2_key 为默认值，跳过合并列。")
  }
} else {
  message("meta.data 中缺少指定的列，跳过合并列。")
}
dev.off()
#save dealed rds
saveRDS(obj, file = out_rds)