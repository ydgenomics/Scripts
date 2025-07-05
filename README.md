# Scripts
*First important work is suming some pipelines and get a gold standard!*
![ydgenomics logo](pngs/ydgenomics_bar_logo.png)

---

## Undering building +++++

### [scATAC](scATAC): *Learning new tool about handing scATAC-seq data*

### [comparsion-RNAVelo](Trajectory/Velocity)

### [Velocity and Trajectory](Trajectory): *Two methods to reveal developmental process*
  Reference: [单细胞RNA速率极简教程 (3)](https://mp.weixin.qq.com/s/JAVNLCZGJlmDkzHwoD106g)

---

## Built
  > If you interested in various GRNs, you could see a a repo named [GRNs](https://github.com/ydgenomics/GRNs)
  > ### [GRN/GRN-IReNA](GRN/GRN-IReNA): *Learning IReNA do combinational analysis of both scRNA and bulk-ATAC*
  > ### [GRN-hdWGCNA](GRN/GRN-hdWGCNA): *Do moduling data and extract gene regulatory network*

  ### [**dataget_scRNAseq**](dataget_scRNAseq): *Perform quality control of single-cell RNA-seq data using SoupX and Scrublet.*

  ### [**DEA-Seurat**](DEA-Seurat): *Conduct differential expression analysis with three Seurat functions to identify gene sets.*

  ### [**enrich_scRNAseq**](enrich_scRNAseq): *Analyze the biological functions of identified genes (e.g., GO and KEGG enrichment).*

  ### [**Ortho_gene**](Ortho_gene): *Compare datasets between two groups to identify orthologous relationships (e.g., using BLASTP).*

  ### [**GRN/GRN-pySCENIC**](GRN/GRN-pySCENIC): *Construct gene regulatory networks with pySCENIC.*

  ### [**cell_similarity**](cell_similarity): *Calculate and visualize similarity across multiple groups using Jaccard, hclust, and MetaNeighbor.*

  ### [**plant_database**](plant_database.md): *Sum plant's database and assess high quality and using databses.*

  ### [**multi_annotation_scRNAseq**](multi_annotation_scRNAseq): *Using various methods of  scRNAseq-data annotation*
---

## Q&A
  - 论先merge再标准化和先标准化再merge的差异: 从可视化的umap来看基本没差异 [codes](temp/test_merge.R)