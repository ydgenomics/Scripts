环境搭建
```shell
conda create -n r r-base=4.2 -y
conda activate r
yum install libicu libicu-devel
conda install conda-forge::r-biocmanager -y
conda install bioconda::bioconductor-clusterprofiler -y
# BiocManager::install("clusterProfiler")
conda install conda-forge::r-tidyverse -y
conda install bioconda::bioconductor-annotationforge -y
```