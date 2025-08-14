# Cluster/Cell similarity with `jaccard`, `hclust` and [metaNeighbor](https://github.com/maggiecrow/MetaNeighbor)
- **Brief:** 不同分组下各个群的相似性
- **Fature**
- **Log:**
  - 250815 增添了`only_metaNeighbor`的判断默认只运行`metaNeighbor`，检查metaNeighbor运行的矩阵对象为`counts`，解决了可视化热图颜色块不显示的问题
- **Tradition:**

---
# Input
- **Variable**
  - `rds` 包含分组信息和cluster信息的rds对象，只跑metaNeighbor只需要rds包含`RNA@counts`,提取`counts`作为singlecellexperiment对象
  - `prefix` 每个`rds`运行之后输出文件的后缀，顺序对应
  - `sankey_order` 每个`rds`里面`batch_key`用','连接用于设定sankey plot的顺序
  - `batch_key` 分组信息键，默认每个`rds`文件的`batch_key`一致
  - `cluster_key` 分群信息键，默认每个`rds`文件的`cluster_key`一致
  - `mem_similarity` 运行similarity的资源(GB)

---
# Output
主要看metaNeighbor结果，hclust和jaccard的结果仅供参考
- **Interpretation**

---
# Detail
- **Overview**

- **Software**
  - **MetaNeighbor**: The output is an AUROC matrix, where each value represents the similarity between cell types across batches or datasets. Higher AUROC values indicate greater similarity and more consistent annotation between datasets.

- **Image**

`as.SingleCellExperiment()` 是 Seurat 提供的方法，用于将 Seurat 对象转换为 SingleCellExperiment 对象（Bioconductor 生态常用的单细胞数据结构）。
转换后的 SingleCellExperiment 对象会保留：
  表达矩阵（默认是 RNA assay 的 data 槽，即 log-normalized 数据）
  细胞元数据（colData，对应 Seurat 的 meta.data）
  基因元数据（rowData，对应 Seurat 的 feature metadata）
  降维结果（如 PCA、UMAP，存储在 reducedDims 中）
```R
# 方法 1：直接提取 counts 并构建 SCE
sce <- SingleCellExperiment(
  assays = list(counts = GetAssayData(sdata, slot = "counts")),
  colData = sdata@meta.data
)

# 方法 2：使用 Seurat::as.SingleCellExperiment() 并指定 slot
sce <- as.SingleCellExperiment(sdata, assay = "RNA", slot = "counts")
```

```R
#使用MetaNeighbor计算每个批次中细胞类型之间的相关性
Aurocs_matrix = MetaNeighborUS(var_genes = global_hvgs, 
                               dat = cca.results.sce, 
                               study_id = cca.results.sce$batch, 
                               cell_type = cca.results.sce$celltype, 
                               fast_version = T)
```

---
# Reference & Citation
> [[R包] MetaNeighbor 第一期 评估不同数据集中细胞类型注释的一致性](https://mp.weixin.qq.com/s/cb9DWJm8zNc1J9wEUNTUVg)
> [常被提起的Jaccard指数是什么？怎么在单细胞中运用和实现Jaccard相似性比较？](https://mp.weixin.qq.com/s/-6iM2phNUh2Qo0wbN0Azpw)
> [练习R：hclust()函数层次聚类分析](https://mp.weixin.qq.com/s/-AvRPX7DG5fzyAVmv8wg7Q)


---
# Coder
- **Editor:** yangdong (yangdong@genomics.cn)
- **GitHub:** [ydgenomics](https://github.com/ydgenomics)
- **Prospect:** Focused on innovative, competitive, open-source projects and collaboration
- **Repository:** [Scripts/enrich_scRNAseq](https://github.com/ydgenomics/Scripts/tree/main/enrich_scRNAseq)
---
