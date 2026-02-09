# 一、单细胞转录组测序技术及其不足
编辑：杨东（yangdong@genomics.cn）
日期：2026/01/15-

---

## 目录
1. 什么是转录组及其测序技术
2. 转录组测序技术的发展
3. 单细胞转录组测序技术
4. 单细胞转录组测序的不足

---

## 1.什么是转录组及其测序技术

- 转录组：对组织/个体所有转录本(各种RNAs)的定性定量分析，因其测序技术的发展，其高通量的优势使其产生了大量的转录组数据，如何解读这些数据以及数据间建立联系值得我们思考和关注。
- 测序技术：转录组测序（FL-cDNA; EST; SAGE, MPSS, CAGE, PET; RNA-Seq）；转录本杂交技术（cDNA-array; Oligonucleotide array; Genome tiling array）

---

## 2.转录组测序技术的发展

- 历史发展
- 技术对比

---

## 3.单细胞转录组测序技术

- 历史发展
- 应用领域
- 技术优势

---

## 4.单细胞转录组测序的不足

- 技术难点
- 数据特征
- 分析困难

# 二、从测序数据到定量矩阵

# 三、定量矩阵去除环境污染RNA

编辑：杨东（yangdong@genomics.cn）
日期：2026/01/15-

---
## 目录

- 环境污染RNA产生的原因
- 现有技术

---

## CAR: Correction of Ambient RNA

```shell
source /opt/software/miniforge/bin/activate
mamba create -n CAR r-base=4.4 -c conda-forge -y
# mamba search r-base
# install.packages('BiocManager')
# BiocManager::install("anndataR")
conda activate CAR
conda install conda-forge::r-seurat -y
conda install bioconda::bioconductor-decontx -y
conda install conda-forge::r-soupx -y

```

---

## Ref:
- 2020_Decontamination of ambient RNA in single-cell RNA-seq with DecontX
- 2020_SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data
- 2023_FastCAR-fast correction for ambient RNA to facilitate differential gene expression analysis in single-cell RNA-sequencing datasets
- 2023_The effect of background noise and its removal on the analysis of single-cell expression data
- 2023_Unsupervised removal of systematic background noise from droplet-based single-cell experiments using CellBender
- 2024_scCDC a computational method for gene-specific contamination detection and correction in single-cell and single-nucleus RNA-seq data
- 2025_Mitigating ambient RNA and doublets effects on single cell transcriptomics analysis in cancer research

---

# 四、去除双细胞

---

## FD: Filter Doublet

---

# 五、细胞分群

---

## Cluster

---

# 六、细胞注释：从细胞群到细胞类型

---

## CA: Cell Annotation

---

# 单细胞数据转换

---

## anndataR
- anndataR | 一行代码读取h5ad为seurat对象 https://mp.weixin.qq.com/s/vN4YxPyF1I0-qTRRD4m16w https://anndatar.data-intuitive.com/index.html

---

# 系统参考文献

- https://www.sc-best-practices.org/preamble.html
- scPlant：一款分析植物单细胞转录组数据的通用工具 https://mp.weixin.qq.com/s/ZEm84pn_3YD7s2CpP3fbpw