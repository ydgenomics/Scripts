- 植物单细胞marker系列
    - [*联川生物*·植物细胞marker数据库总览，植物单细胞分析的最佳伴侣！| 植物单细胞专题](https://mp.weixin.qq.com/s/CXGkNuBDQin5MrPWMgt8ng)
    - [scPlantDB](https://biobigdata.nju.edu.cn/scplantdb/home) [*基迪奥生物*·分享一个好用的植物单细胞数据库](https://mp.weixin.qq.com/s/1dTCDc5U3dvCy15GfLRY4A)
    - [PlantCellMarker](https://www.tobaccodb.org/pcmdb/homePage) [*生信益站*·单细胞专题25| 植物细胞类型注释数据库: PlantCellMarker](https://mp.weixin.qq.com/s/Y1AyXa8jkQBV4yWo_HihTw)
    - [PsctH](http://jinlab.hzau.edu.cn/PsctH/) [*植物科学最前言*·PBJ | 华中农大开发出植物单细胞转录组综合数据库，提供综合全面的单细胞Marker基因资源和单细胞研究的workflow](https://mp.weixin.qq.com/s/5dMORWQeX4eTFgH0e1YkTg)
    - [PlantscRNAdb](http://ibi.zju.edu.cn/plantscrnadb/index.php)

第一种策略：利用现有的基因集和参考数据库
第二种策略：借助人工智能（AI）工具 scGPT scPlantLLM
第三种策略：应用机器学习方法
第四种策略：使用算法方法或专用软件
第五种策略：结合领域知识进行手动注释


测试
![scPlantDB下载拟南芥根的1w细胞数据集做测试](png/download_testdata.png)
```R
# 无法顺利构建对象，那我就先转h5ad来处理这个问题
seu <- readRDS("/data/work/multi_anno/AT_root_SRP273996.rds")
seu
# An object of class Seurat 
# 79128 features across 10431 samples within 3 assays 
# Active assay: integrated (2000 features, 2000 variable features)
#  2 layers present: data, scale.data
#  2 other assays present: RNA, SCT
#  2 dimensional reductions calculated: pca, umap
DefaultAssay(seu) <- "RNA"; seu
# An object of class Seurat 
# 79128 features across 10431 samples within 3 assays 
# Active assay: RNA (53678 features, 0 variable features)
#  2 layers present: counts, data
#  2 other assays present: SCT, integrated
#  2 dimensional reductions calculated: pca, umap
rna_data <- GetAssayData(seu, assay = "RNA", layer = "counts")
seu2 <- CreateAssayObject(counts = rna_data, name = "RNA")
seu2 <- AddMetaData(seu2, metadata = seu@RNA@meta.data); seu2
```
```R
sceasy::convertFormat(temp0,from = "seurat",to = "anndata",assay = "RNA",main_layer = "counts",outFile = "/data/work/multi_anno/AT_root_SRP273996_RNA_rh.h5ad")
sceasy::convertFormat(temp0,from = "seurat",to = "anndata",assay = "SCT",main_layer = "counts",outFile = "/data/work/multi_anno/AT_root_SRP273996_SCT_rh.h5ad")
sceasy::convertFormat(temp0,from = "seurat",to = "anndata",assay = "integrated",main_layer = "counts",outFile = "/data/work/multi_anno/AT_root_SRP273996_integrated_rh.h5ad")
```
[单细胞全自动注释篇(四)——ScType](https://mp.weixin.qq.com/s/hKBiZCHwDdoJOk0YChbtMA)
sctype：csv格式的marker基因列表。关注查询数据集的scale.data的矩阵，按查询分群来做注释，一个群可能会被分到多个细胞类型，取最优，同时如果太差会被认定为unknown。csv的marker基因要尽可能多的存在于查询数据的基因中

[使用singleR基于自建数据库来自动化注释单细胞转录组亚群](https://mp.weixin.qq.com/s/GpOxe4WLIrBOjbdH5gfyOQ)
singleR: 拿到参考数据集的RNA@counts矩阵后，计算每种细胞的平均表达量后做log处理。而查询数据集的data矩阵按每个细胞去拟合参考数据集的细胞类型表达模式，所有两个数据集共同所有的基因数很重要
我有一个疑问,构建reference对counts做的logcounts处理和查询数据做了Normalization的data，在计算方式上不一样，这样会影响singleR注释。假设不考虑两次数据实验误差。我想到了，如果singleR设计的时候查看的是不同基因的波动而非专注于某个值的话，就不会在意这个问题
logcounts()：对原始计数数据进行对数转换，减少数据的偏态分布。输出对数转换后的矩阵，值的范围通常在0到10之间。
Normalization()：对原始计数数据进行归一化处理，调整每个细胞的总读数。输出归一化后的矩阵，值的范围通常在0到1之间。

# Interpretation of Results
根圆（ord == 1 的 cluster 节点） 的半径 ∝ nodes_lvl1$ncells（即该 cluster 的细胞总数）。
子圆（ord == 2 的 cell-type 节点） 的半径 ∝ nodes_lvl2$ncells（即该类型在该 cluster 中的 scores，你把它赋给了 ncells）