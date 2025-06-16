环境搭建
[clusterprofiler安装-六种方法](https://mp.weixin.qq.com/s/BYydet8hoBdbfZQgdc2dIA)
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
> [基于AnnotationDbi非模式物种OrgDB下载制作](https://mp.weixin.qq.com/s/auyTKJhfos0wi_yPsA7O0g)
> [从 gtf 文件构建 orgdb 和 txdb 数据库](https://mp.weixin.qq.com/s/w3FFimm-xF2OY20aoFRcSg)

**解决了之前要安装包才能调用库的问题**
```R 
orgdb <- loadDb("/data/work/0.peanut/orgdb/output/org.Ahypogaea.eg.db/inst/extdata/org.Ahypogaea.eg.sqlite") #加载本地数据库
keytypes(orgdb)  # 查看这个数据库中有哪几种keytypes
#  [1] "EVIDENCE"    "EVIDENCEALL" "GENENAME"    "GID"         "GO"         
#  [6] "GOALL"       "Ko"          "ONTOLOGY"    "ONTOLOGYALL" "Pathway"    
length(keys(orgdb)) #查看包含的基因数量
# [1] 68781
columns(orgdb)  #查看OrgDb对象的数据类型
#  [1] "EVIDENCE"    "EVIDENCEALL" "GENENAME"    "GID"         "GO"         
#  [6] "GOALL"       "Ko"          "ONTOLOGY"    "ONTOLOGYALL" "Pathway" 
saveDb(orgdb,file="/data/work/0.peanut/orgdb/output/Ahypogaea.Orgdb") #把Capra_hircus对象保存成Capra_hircus.OrgDb文件。
```

enrichplot可视化和gsea分析后续补上