### [SCPiplines](): *Pack some codes as functions to high efficiently use*
  Although SCP had do this, but it could do personal works. If I want do some new works, it couldn't adapt to my need, so i decided do myself single cell piplines.
  - First question is building a image of enough do most works
```shell
# R of environment: SCr
conda create -n CHOIR r-base=4.2.2 -y
conda activate CHOIR
# devtools::install_github("pengminshi/mrtree")
#conda install icu=58.2
#devtools::install_github("pengminshi/mrtree")

conda create -n SCr r-base=4.3 -y
conda activate SCr
conda install conda-forge::r-biocmanager -y
conda install conda-forge::r-devtools -y
conda install conda-forge::r-remotes
conda install conda-forge::r-seurat -y
#Rscript -e 'remotes::install_github("corceslab/CHOIR", ref="main", repos = BiocManager::repositories(), upgrade = "never")'
Rscript -e 'remotes::install_github("corceslab/CHOIR", ref="dev", repos = BiocManager::repositories(), upgrade = "never")'
conda install conda-forge::r-soupx -y
conda install bioconda::bioconductor-decontx -y
conda install conda-forge::r-hgnchelper -y
conda install bioconda::bioconductor-singler -y
conda install conda-forge::r-harmony -y
#conda install bioconda::r-presto -y
# library(devtools)
# install_github('immunogenomics/presto')
conda install bioconda::bioconductor-clusterprofiler -y
conda install bioconda::bioconductor-aucell -y
Rscript -e 'install.packages("GeneNMF")'
Rscript -e 'install_github("Jasonxu0109/PlantPhoneDB")'
Rscript -e 'BiocManager::install(c("WGCNA", "UCell", "GenomicRanges", "GeneOverlap"))'
conda install bioconda::bioconductor-metaneighbor -y
conda install conda-forge::r-optparse -y
conda install bioconda::bioconductor-dropletutils -y
conda install conda-forge::r-tidyverse -y
conda install bioconda::bioconductor-annotationforge -y
# R plot
conda install conda-forge::r-openxlsx -y
conda install bioconda::bioconductor-complexheatmap -y
conda install r::r-ambient -y
conda install conda-forge::r-ggsci -y
```
```R
library(Seurat)
library(CHOIR)
library(presto)
library(SoupX)
library(decontX)
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
library(SingleR)
library(harmony)
library(clusterProfiler)
library(AUCell)
library(GeneNMF)
library(PlantPhoneDB)
library(WGCNA)
library(ComplexHeatmap)
library(MetaNeighbor)
library(optparse)
library(DropletUtils)
library(tidyverse)
library(AnnotationForge)
```

```shell
# python of environment: SCpy
conda create -n SCpy python=3.12 -y
conda install -c conda-forge scanpy python-igraph leidenalg -y
pip install memento-de
conda install bioconda::scrublet -y
```