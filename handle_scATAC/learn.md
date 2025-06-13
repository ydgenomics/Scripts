cotton项目
数据有bulk RNA； bulk ATAC； scRNAseq的数据。单细胞技术发展为什么还要做bulk sequence？因为single cell水平的转录本捕获低，或者说基因面覆盖的广，bulk保证了其深度！

这个项目是基于时序数据推断调控网络，从分子的调控网络去解释生物的生长发育的事件。
如何去联合分析？

IReNA环境配置
```shell
conda create -n IReNA r-base=4.0 -y
conda activate IReNA
conda install bioconda::bioconductor-rsamtools -y
conda install bioconda::bioconductor-chipseeker -y
conda install bioconda::bioconductor-monocle -y
conda install bioconda::bioconductor-rcistarget -y
conda install bioconda::bioconductor-rcy3 -y
conda install bioconda::bioconductor-clusterprofiler -y
conda install conda-forge::r-devtools -y
```
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.12")

BiocManager::install(c('Rsamtools', 'ChIPseeker', 'monocle',
                       'RcisTarget', 'RCy3', 'clusterProfiler'))

install.packages("devtools")
devtools::install_github("jiang-junyao/IReNA")

library(IReNA)
```

[pp.bulkATAC](https://jiang-junyao.github.io/IReNA/ATAC-seq-preprocessing)
[pp.scRNAseq](https://jiang-junyao.github.io/IReNA/scRNA-seq-preprocessing)
[bulkATAC+scRNAseq](https://jiang-junyao.github.io/IReNA/scATAC+scRNA)