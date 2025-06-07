# Date: 20250607
# Reference: [*生信钱同学*·全代码干货奉上——多样本多方案去除单细胞环境RNA污染——这次把这个聊清楚](https://mp.weixin.qq.com/s/1eJq3u-aKpQaL9CM7bV94g)

samples <- dir(path = "path_to_cellranger")
cl <- makeCluster(0.75 * detectCores())

data.tod <- parLapply(
  cl, samples,
  function(sample) {
	if (!require(Seurat)) install.packages("Seurat")
	path <- paste0(
	  "path_to_cellranger",
	  paste0(sample, "/raw_feature_bc_matrix")
	)
	CreateSeuratObject(
	  counts = Read10X(data.dir = path),
	  project = sample,
	  min.cells = 0,
	  min.features = 0
	)
  }
)
names(data.tod) <- samples

data.toc <- parLapply(
  cl, samples,
  function(sample) {
	if (!require(Seurat)) install.packages("Seurat")
	path <- paste0(
	  "path_to_cellranger",
	  paste0(sample, "/filtered_feature_bc_matrix")
	)
	CreateSeuratObject(
	  counts = Read10X(data.dir = path),
	  project = sample
	)
  }
)
names(data.toc) <- samples

preprocess <- function(sce) {
  sce <- CreateSeuratObject(sce)
  sce <- NormalizeData(sce)
  sce <- FindVariableFeatures(sce)
  sce <- ScaleData(sce)
  sce <- RunPCA(sce)
  sce <- FindNeighbors(sce, dims = 1:20)
  sce <- FindClusters(sce, resolution = 0.5)
  sce <- RunUMAP(sce, dims = 1:20)
  return(sce)
}

pp.toc <- lapply(
  1:length(data.toc),
  function(i) {
	sce = data.toc[[i]]
	sce = preprocess(sce)
	return(sce)
  }
)

######### SoupX #########
library(SoupX)
library(Seurat)
library(Matrix)
# 定义SoupX运行函数，输出矫正矩阵
run_soupx <- function(toc, tod) {
  sc = SoupChannel(tod, toc)
  sc = setClusters(sc, setNames(sce$seurat_clusters, colnames(sce)))
  sc = autoEstCont(sc)
  out = adjustCounts(sc)
  return(out)
}
# 从构建好的seurat list中提取过滤矩阵toc和原始矩阵tod
out.list = lapply(
  1:length(pp.toc),
  function(i) {
    toc = GetAssayData(pp.toc[[i]])
    tod = GetAssayData(data.tod[[i]])
    temp = run_soupx(toc, tod)
    return(temp)
  }
)
# 合并矫正矩阵，转换为seurat对象
out = do.call(cbind, out.list)
sce = CreateSeuratObject(count = out)
rm(out.list)
rm(out)

####### decontx #########
# decontX运行需要基于sce对象进行，首先将构建好的seurat list合并为完整的
# seurat，然后转换为sce对象。
filter <- merge(
  data.toc[[1]],
  y = data.toc[-1],
  add.cell.ids = samples
)
sce.filter <- as.SingleCellExperiment(filter)

raw <- merge(
  data.tod[[1]],
  y = data.tod[-1],
  add.cell.ids = samples
)
sce.raw <- as.SingleCellExperiment(raw)

# 添加batch参数指定批次和原始背景计数（raw）
sce <- decontX(
  sce.filter,
  batch = sce.filter$orig.ident,
  background = sce.raw,
  bgBatch = sce.raw$orig.ident
)

# 提取污染分数和矫正矩阵
filter$contamination <- sce$decontX_contamination
filter[["decontX"]] <- CreateAssay5Object(counts = assays(sce)[[2]])
sce <- filter

# 没有空液滴过滤（适合一些只给了不跟数据的公共数据挖掘）
sce.filter <- as.SingleCellExperiment(filter)
sce <- decontX(sce.filter, batch = sce.filter$orig.ident)
filter$contamination <- sce$decontX_contamination
filter[["decontX"]] <- CreateAssay5Object(counts = assays(sce)[[2]])
sce <- filter


####### scCDC[3]（刚发的一篇文章，也可以试试） #####
