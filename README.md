# Scripts
*First important work is sum all workflow and get a gold standard!*

## Undering building……
## [**plant_database**](): *Sum plant's database and assess high quality and using databses.*
  - Reference:
    - [*小杜的生信笔记*·植物学中常用的数据库 | 通用数据库](https://mp.weixin.qq.com/s/eWRKpZbVN8iY1qmu5mue2g)
    - [*基迪奥生物*·研究植物转录调控，你不能不知道的数据库](https://mp.weixin.qq.com/s/yee680uNUmQQUOXISr479A) [PlantTFDB](http://planttfdb.cbi.pku.edu.cn/)
    - 植物单细胞marker系列
      - [*联川生物*·植物细胞marker数据库总览，植物单细胞分析的最佳伴侣！| 植物单细胞专题](https://mp.weixin.qq.com/s/CXGkNuBDQin5MrPWMgt8ng)
      - [scPlantDB](https://biobigdata.nju.edu.cn/scplantdb/home) [*基迪奥生物*·分享一个好用的植物单细胞数据库](https://mp.weixin.qq.com/s/1dTCDc5U3dvCy15GfLRY4A)
      - [PlantCellMarker](https://www.tobaccodb.org/pcmdb/homePage) [*生信益站*·单细胞专题25| 植物细胞类型注释数据库: PlantCellMarker](https://mp.weixin.qq.com/s/Y1AyXa8jkQBV4yWo_HihTw)
      - [PsctH](http://jinlab.hzau.edu.cn/PsctH/) [*植物科学最前言*·PBJ | 华中农大开发出植物单细胞转录组综合数据库，提供综合全面的单细胞Marker基因资源和单细胞研究的workflow](https://mp.weixin.qq.com/s/5dMORWQeX4eTFgH0e1YkTg)
      - [PlantscRNAdb](http://ibi.zju.edu.cn/plantscrnadb/index.php)

## Built
  ### [**dataget_scRNAseq**](dataget_scRNAseq): *Perform quality control of single-cell RNA-seq data using SoupX and Scrublet.*

  ### [**DEA-Seurat**](DEA-Seurat): *Conduct differential expression analysis with three Seurat functions to identify gene sets.*

  ### [**enrich_scRNAseq**](enrich_scRNAseq): *Analyze the biological functions of identified genes (e.g., GO and KEGG enrichment).*

  ### [**Ortho_gene**](Ortho_gene): *Compare datasets between two groups to identify orthologous relationships (e.g., using BLASTP).*

  ### [**GRN-pySCENIC**](GRN-pySCENIC): *Construct gene regulatory networks with pySCENIC.*

  ### [**cell_similarity**](cell_similarity): *Calculate and visualize similarity across multiple groups using Jaccard, hclust, and MetaNeighbor.*