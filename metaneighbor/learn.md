# Learn metaneighbor
[github metaneighbor](https://github.com/maggiecrow/MetaNeighbor)
[[R包] MetaNeighbor 第一期 评估不同数据集中细胞类型注释的一致性](https://mp.weixin.qq.com/s/cb9DWJm8zNc1J9wEUNTUVg)
包含应用的两个场景和结果解读

场景一：与已知注释的数据集进行相关性分析这种方法是通过比较目标数据集与已知、良好注释的数据集进行相关性分析，来验证注释的准确性。例如，在研究非模式植物时，可以采取以下步骤：基因比对：将非模式植物的基因与已知模式植物（如拟南芥）的基因进行比对。使用已知标记基因进行注释：利用模式植物中已验证的细胞类型标记基因对自己的数据进行初步注释。相关性分析：使用表达数据分析非模式植物与模式植物之间相同细胞类型的表达模式。注释到的相同细胞类型之间的表达相关性越高，表明注释的准确性越理想。
场景二：使用特定标记基因进行人工注释当缺乏已知的参考数据集时，可以使用特定于研究对象的标记基因进行人工注释：标记基因选择：选择非模式植物特有的、已知功能的标记基因。人工注释：使用这些标记基因对数据进行人工注释，确定细胞类型。跨数据集验证：比如我研究一种新的植物-菜花！现有两个菜花数据集，我在这两个数据集中都能注释出相同的细胞亚型（例如，细胞类型A、B、C）。在各个数据集中注释出的相同标签的细胞类型应表现出高度的表达相似性，这有助于证明注释的稳健性和一致性。
这里计算的是矩阵，应该也涉及到reference数据基因名和test数据基因名覆盖度的问题

```R
#使用MetaNeighbor计算每个批次中细胞类型之间的相关性
Aurocs_matrix = MetaNeighborUS(var_genes = global_hvgs, 
                               dat = cca.results.sce, 
                               study_id = cca.results.sce$batch, 
                               cell_type = cca.results.sce$celltype, 
                               fast_version = T)
```
计算的是矩阵的那一层呢？counts还是data
