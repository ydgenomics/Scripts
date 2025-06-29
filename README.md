# Scripts
*First important work is sum all workflow and get a gold standard!*

## Q&A
  - 论先merge再标准化和先标准化再merge的差异

## Undering building……
> If you interested in various GRNs, you could see a a repo named [GRNs](https://github.com/ydgenomics/GRNs)

### [GRN/GRN-IReNA](GRN/GRN-IReNA): *Learning IReNA do combinational analysis of both scRNA and bulk-ATAC*
### [GRN-hdWGCNA](GRN/GRN-hdWGCNA): *Do moduling data and extract gene regulatory network*

### [scATAC](scATAC): *Learning new tool about handing scATAC-seq data*

### [comparsion-RNAVelo](Trajectory/Velocity)

### [Velocity and Trajectory](Trajectory): *Two methods to reveal developmental process*
  Reference: [单细胞RNA速率极简教程 (3)](https://mp.weixin.qq.com/s/JAVNLCZGJlmDkzHwoD106g)

### [test_data](): *~~Test data of model organism for check all analysis pipline~~*
  ~~[The data of *Arabidopsis thaliana*(True leaf)](test_data/ERP132245.h5ad) downloaded from scPlantdb, having 2018 cells and two conditions including mild drought and normal.~~

### [Project help learn: *Peanut* and *AT*]: I wish use the high quality data as test, building my workflow and improve my understanding for codes.
  - 0608 change its gene names of matrix
  - 0607 run dataget_scRNAseq

## Built
  ### [**dataget_scRNAseq**](dataget_scRNAseq): *Perform quality control of single-cell RNA-seq data using SoupX and Scrublet.*

  ### [**DEA-Seurat**](DEA-Seurat): *Conduct differential expression analysis with three Seurat functions to identify gene sets.*

  ### [**enrich_scRNAseq**](enrich_scRNAseq): *Analyze the biological functions of identified genes (e.g., GO and KEGG enrichment).*

  ### [**Ortho_gene**](Ortho_gene): *Compare datasets between two groups to identify orthologous relationships (e.g., using BLASTP).*

  ### [**GRN/GRN-pySCENIC**](GRN/GRN-pySCENIC): *Construct gene regulatory networks with pySCENIC.*

  ### [**cell_similarity**](cell_similarity): *Calculate and visualize similarity across multiple groups using Jaccard, hclust, and MetaNeighbor.*

  ### [**plant_database**](plant_database.md): *Sum plant's database and assess high quality and using databses.*

  ### [**multi_annotation_scRNAseq**](multi_annotation_scRNAseq): *Using various methods of  scRNAseq-data annotation*

