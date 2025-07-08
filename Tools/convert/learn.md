# 精简convert代码，让其更容易调用，更通用

```shell
source /opt/software/miniconda3/bin/activate
conda create -n sc r-base=4.3 -y
conda activate sc
conda install conda-forge::r-seurat -y
conda install conda-forge::r-reticulate -y
conda install bioconda::r-sceasy -y
conda install conda-forge::r-devtools -y
# devtools::install_github("cellgeni/schard")
conda install conda-forge::scanpy -y
conda install conda-forge::loompy -y
```