[ATAC中文学习资源](https://mp.weixin.qq.com/mp/appmsgalbum?__biz=MzAxMDkxODM1Ng==&action=getalbum&album_id=3825619502127398912&subscene=189&scenenote=https%3A%2F%2Fmp.weixin.qq.com%2Fs%3F__biz%3DMzAxMDkxODM1Ng%3D%3D%26mid%3D2247538904%26idx%3D1%26sn%3Db8a8393137c5e599e5f7a1b33715d4e9%26chksm%3D9b4b1863ac3c91750db3d0d5f0732905197bffcf786454309094cebc96e797c07f6a04dd1fb0%26cur_album_id%3D3825619502127398912%26scene%3D189%26key%3Ddaf9bdc5abc4e8d0851611211b85885dbdc306f6d992b3cef8eef0c0bcbdab0f778f294c8cf6dc082f581bcc7372b6ba2704cfeba7cc937b4c98742ff897e7d4ef4aa9ebe4ec334af028e6b504b6ebaf6884621c2102adb09d725771c46bffe01085ee6b39742271db5cf8ed655d013022afe79bf821c1cab206482948f7071b%26ascene%3D0%26uin%3DNDIxMzk4MTk3%26devicetype%3DWindows%2B11%2Bx64%26version%3D63090c33%26lang%3Dzh_CN%26countrycode%3DCN%26exportkey%3Dn_ChQIAhIQj6Pe6cyNMp0D49NAt6LmZhLmAQIE97dBBAEAAAAAADdcNmS%252Ba1sAAAAOpnltbLcz9gKNyK89dVj0i0JJwIhsKgwUiurYkjY88FVyEcYdGphnQEO%252BTD7vXb3GyCJ08cCSqMxRu84V%252BKSftWd8WdCpmr2pJkv%252BhE2vorTlAHRfqzfloGOnFLfdITdFvP1MHIjPdUzj%252BVCQuDfaSkJ0KH4hjvpgjupLYHVsnw1uIufC9PqVNG2pMTYP4hofIoopMBQ9F9hhjKYkscLj9hBb2Ax8mwR1jr1PmdfM37AmfZb01SkUYHoKOpeyRqTToCjQBequnmoNVvx9fTxj%26acctmode%3D0%26pass_ticket%3D%252BzkNfT%252FzC8z%252F2ygWf1AwYZ%252FAlE81E1r%252BDRKfc8SQYA1gtHVgExFF1b108S4Q5jGj%26wx_header%3D1&nolastread=1&sessionid=1742267548#wechat_redirect)

[1.综述：ATAC-Seq 数据分析工具大全](https://mp.weixin.qq.com/s/DSm7z9Z2QYfZJ5A_Ch7Z6g)

cotton项目
数据有bulk RNA； bulk ATAC； scRNAseq的数据。单细胞技术发展为什么还要做bulk sequence？因为single cell水平的转录本捕获低，或者说基因面覆盖的广，bulk保证了其深度！

这个项目是基于时序数据推断调控网络，从分子的调控网络去解释生物的生长发育的事件。
如何去联合分析？

> IReNA环境配置
```shell
conda create -n IReNA r-base=4.3 -y
conda activate IReNA
conda install bioconda::bioconductor-rsamtools -y
conda install bioconda::bioconductor-chipseeker -y
conda install bioconda::bioconductor-monocle -y
conda install bioconda::bioconductor-rcistarget -y
conda install bioconda::bioconductor-rcy3 -y
conda install bioconda::bioconductor-clusterprofiler -y
conda install conda-forge::r-devtools -y
conda install conda-forge::r-remotes -y
conda install conda-forge::r-pbapply -y
conda install conda-forge::r-rocr -y
conda install conda-forge::r-seurat -y # 5.3
conda install bioconda::bioconductor-edger -y
conda install bioconda::bioconductor-edger -y
conda install conda-forge::r-furrr -y
conda install conda-forge::r-future -y
```

```R
devtools::install_github("jiang-junyao/IReNA")
library(IReNA)
library(pheatmap)
library(DDRTree)
library(RcisTarget)
library(Rsamtools)
library(GenomicRanges)
library(ChIPseeker)
library(IRanges)
library(Biostrings)
library(igraph)
library(pbapply) # conda install conda-forge::r-pbapply -y
library(RCy3)
library(ROCR) # conda install conda-forge::r-rocr -y
library(Seurat) # conda install conda-forge::r-seurat -y # 5.3
library(ggplot2)
library(dplyr)
library(gridExtra)
library(stats)
library(utils)
library(monocle)
library(VGAM)
library(BiocGenerics)
library(clusterProfiler)
library(edgeR) # conda install bioconda::bioconductor-edger -y
library(limma)
library(rlang)
library(reshape2)
library(stringr)
library(furrr) # conda install conda-forge::r-furrr -y
library(future) # conda install conda-forge::r-future -y
library(purrr)
```

[pp.bulkATAC](https://jiang-junyao.github.io/IReNA/ATAC-seq-preprocessing)
```shell
conda create -n macs3
conda activate macs3
conda install -c bioconda macs3 -y
```

```shell
conda install bioconda::htseq -y
```
[pp.scRNAseq](https://jiang-junyao.github.io/IReNA/scRNA-seq-preprocessing)

## [bulkATAC+scRNAseq](https://jiang-junyao.github.io/IReNA/scATAC+scRNA)
[samtools]() *Image: chromap*
[MACS3](https://macs3-project.github.io/MACS/index.html) *Image: macs3*
[IReNA]() *Image: IReNA *
[htseq](https://github.com/htseq/htseq) *Image: htseq*
[hint](https://reg-gen.readthedocs.io/en/latest/hint/introduction.html) *Image: rgt-hint*
