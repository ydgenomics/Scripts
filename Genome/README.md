# The Annotation of Genome

---
数据质控-组装-注释
**THINK:** 基因组组装和注释，测序拿到的序列为*read*，*read*分长短，二代为短，三代测序为长，通过加测序深度可以找到*read*间重叠部分，进而实现*read*间的连接得到整个基因序列。在做这一步之前应该还有*read*质控，并不是所有的*read*都会在这个组装考虑之中，所以什么样的的*read*才算可以作为考虑的*read* [测序数据质量控制:一键完成质量过滤、去重复与序列修剪](https://mp.weixin.qq.com/s/qgcYSrUwyOFAlFwfp1_Knw) [fastq和fasta格式文件](https://mp.weixin.qq.com/s/t-Z0S-1ByqAWAf3IZ4YA4w) [二、RNA-seq数据分析系列教程：FASTA与FASTQ文件格式解读学习](https://mp.weixin.qq.com/s/zAjplLsCgvDuydxxkNEksQ) [学习笔记丨fasta与fastq格式相互转换，使用Python脚本轻松实现！](https://mp.weixin.qq.com/s/POE4EkSmnK7fMcNDjpWCcw) 。数据组装分有参和无参两种 [基因组组装全攻略：从线粒体到全基因组](https://mp.weixin.qq.com/s/zmLAgaA96SKVl10xouSnYA)，构建T2T

PROBLEM:
- gtf转gff
- EVM是否可以只支持两种输入的gff


**Q&A:**
  - 测序深度

[生物学功能注释三板斧](https://mp.weixin.qq.com/s/YVNOeYEIvh7adXb4nSOO5A)
[braker3--基因结构注释](https://mp.weixin.qq.com/s/R9kxyTfjkZ_iTs0mqwNFKw)
[一个染色体组级别的基因组注释全流程代码分享](https://mp.weixin.qq.com/s/ZFAbzxRgf4E0OqRQLkte5w)

# Reference & Citation
- 序列名词 [基因序列中的一些名词区别（CDS、Exon、Intron、UTR、ORF、启动子、TF等）](https://zhuanlan.zhihu.com/p/557609219)
- [LiftOn 基因组注释迁移工具](https://mp.weixin.qq.com/s/SZEFfG8Q6qwlj9Fxb_ostw)
- [Liftoff：基于参考基因组的基因组注释](https://mp.weixin.qq.com/s/jiG1QVBQQfhdSMzmucCIJw)
- [工具分享|agat:高效玩转GTF/GFF基因注释文件](https://mp.weixin.qq.com/s/U14Y-qAqvh9BjyiG1nSR1A) 
- [website ATAG](https://agat.readthedocs.io/en/latest/index.html)

- [使用OMArk评估基因组注释](https://mp.weixin.qq.com/s/eacxbwZQhXtZM8BJ6BS3dg)
- [提取最长转录本以及用BUSCO评估基因注释的完整性](https://mp.weixin.qq.com/s/hiOtkY80kQQFKd2GjqbRkQ)

- [从零开始基因组注释（七）| EVidenceModeler 基因结构证据整合](https://mp.weixin.qq.com/s/q-bEHXu8Vu0Y8hUPwuAvrA)
- [EVidenceModeler整合基因组注释结果](https://mp.weixin.qq.com/s/3Brjr8mBnP7O70ba0aABJQ)
- [关于基因结构注释的一些感悟](http://xuzhougeng.com/archives/reflections-on-gene-structure-annotation)