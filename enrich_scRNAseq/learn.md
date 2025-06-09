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

gofigure的背景文件下载`ic.tsv`, `relations_full.tsv`, `go.obo` [gofigure rep](https://gitlab.com/evogenlab/GO-Figure/-/tree/master/data?ref_type=heads) [buidu cloud pwd:1111]()

# Reference
> [模式植物GO背景基因集制作](https://mp.weixin.qq.com/s/08hAZs24mi_KBOa4QZRLdQ)
> [模式植物构建orgDb数据库 | 以org.Slycompersicum.eg.db为例](https://mp.weixin.qq.com/s/b8OrDKJJGdXwF9B1C7l6zg)