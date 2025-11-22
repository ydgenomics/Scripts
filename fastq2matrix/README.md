# fastq2matrix: Embrace shared data

---
关于从fastq到拿到矩阵。fastq文件是什么，与fasta的区别是什么。fastq与fasta格式很像，但是比fasta更复杂一点，会多一行碱基质量分数，怎么来的呢(不知道)。fasta的格式为序列名和核酸序列，将fastq去除质量分数就可以那到fasta。如何从测序read到fasta/fastq，测序拿到的都是一些片段，基于不同的测序技术其片段大小各异，但是一个片段很难是一个基因，要想拿到一个基因的序列，需要将这些片段做组装再注释。这里面就会有基因组的知识。[RNAseq上游分析全流程详解](https://mp.weixin.qq.com/s/ieqYQwabaYzu62wC6z9stQ)
- **Image**
  - cellranger
  - fastq2matrix

```shell
cd /opt
# [ download file from downloads page ]
wget -O cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1757082956&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=EhdqMA~vIQizKQIIkf9QIF04dpc3bWDQZm7ehNWxui~Ne1KOh0Zbs8Vchipe7hHrCD36bSc1Cn5u5E8gkv12ILO~xcZKpLlwOmYYkkz3IJGYPXkaTcLjCVW9LXGs1iAm7cNrGoeuPd~ZV1csplqx~aXytb0Kdc-RAcNWPf~2Gr2uf18yMdEDaZ~zpQq~ZaaG09yirqM-wJ9oyi0B6aNN5xNIg9qsVb63R3vdFcy0i2bMhKH~L9~6Hw9jEyuA4E5jXO98LUCwz2JX7qtN0ioWr9wRPX4GlHqScBh4cz-tWQT3NbBclbK~XRULrSNLdHxA0LQO-xk6ADN4JfCDHiWS4w__"
tar -xzvf cellranger-9.0.1.tar.gz
export PATH=/opt/cellranger-9.0.1:$PATH # /opt/cellranger-9.0.1/bin/cellranger
rm cellranger-9.0.1.tar.gz
source /opt/software/miniconda3/bin/activate
conda install bioconda::sra-tools -y
# export PATH=/opt/software/miniconda3/bin:$PATH # 
conda install conda-forge::pigz -y
```

<details> <summary><strong> 棉花数据复现 </strong></summary>

# [Single-cell transcriptome atlas identified novel regulators for pigment gland morphogenesis in cotton](https://doi.org/10.1111/pbi.14035)
陆地棉子叶（1周）的单细胞测序 。 Glanded代表CCRI12栽培种陆地棉，Glandless代表CCRI12突变体（无腺体）。 10X Genomics
  - SRR31330970 Glanded
  - SRR31330969 Glandless

```shell
cd /data/work/Reference/Gossypium
prefetch SRR31330970 --max-size 120G
id="SRR31330970"
fastq-dump -O fastq/ --split-3 --gzip ./${id}/${id}.sra
fasterq-dump -O fastq/${id} --split-files -e 40 ./${id}/${id}.sra  --include-technical  -x
prefetch SRR31330969 --max-size 120G

# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/990/345/GCF_007990345.1_Gossypium_hirsutum_v2.1/GCF_007990345.1_Gossypium_hirsutum_v2.1_genomic.gtf.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/990/345/GCF_007990345.1_Gossypium_hirsutum_v2.1/GCF_007990345.1_Gossypium_hirsutum_v2.1_protein.faa.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/990/345/GCF_007990345.1_Gossypium_hirsutum_v2.1/GCF_007990345.1_Gossypium_hirsutum_v2.1_genomic.fna.gz
```

复现过程找到原文的参考基因组很困难，通过NCBI的refgenome来找，你们的基因名都是LOC命名方式，这些只是标识符，从标识符到相比更具意义的基因名有不少困难

pigment gland 色素腺体
cottonseed 棉籽
gossypol 棉子酚
derivative 派生物/衍生物
gland cotton ‘CCRI12’ and glandless cotton ‘CCRI12gl’
cotyledon 子叶
protoplast 原生质体

A total of 9186 individual cells, including 4790 cells from ‘CCRI12’ and 4396 cells from ‘CCRI12gl’, were obtained after cell filtering process (Figure S1, Table S1) and were divided into 12 clusters based on highly variable genes (Figure 1b, Figure S2).

  
</details>

---
# Reference & Citation
- 数据下载
  - 技能树在文末罗列了自己写的cellranger[删一半的细胞后降维聚类分群当然好看啊](https://mp.weixin.qq.com/s/xVAkPRq_b1g1Sj7jDseoxQ)
  - [单细胞上游分析（一）：单细胞测序原始数据下载和解压](https://mp.weixin.qq.com/s/49zzxtWwu6RMY7g79dsV6Q)
  - [2. NCBI下载文献序列数据，并基于SRA Toolkit转化数据为fastq文件](https://mp.weixin.qq.com/s/WZ1fr1IB6ngDlaCrP0zt-g)
  - [新手别愁！超详细的NCBI数据下载流程](https://mp.weixin.qq.com/s/CWGdJDUsd1_jRCbnH0DOmw)
  - [https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)
  - [使用pigz快速压缩当前文件夹下面的全部fastq测序数据](https://mp.weixin.qq.com/s/4pWkIoUAKvFoZrYKoDHLJg)
  - [NCBI下载的sra文件的拆分](https://mp.weixin.qq.com/s/rWLr4Xp8weztknznjL8XZQ)
- cellranger
  - [https://www.10xgenomics.com/support/cn/software/cell-ranger/downloads](https://www.10xgenomics.com/support/cn/software/cell-ranger/downloads)
  - [https://www.10xgenomics.com/support/cn/software/cell-ranger-arc/downloads#download-links](https://www.10xgenomics.com/support/cn/software/cell-ranger-arc/downloads#download-links)
  - [单细胞转录组：下载 NCBI 数据 （SRA ）：从fastq 到10x矩阵](https://mp.weixin.qq.com/s/4tqmyMEkCtn0UJwN-l4Swg)
  - [单细胞专题31| Cellranger: 10x 转录组从原始下机数据到表达矩阵](https://mp.weixin.qq.com/s/j0iZFRKaQrVV3Lq4TVzWuQ)
  - [快速上手Cellranger](https://mp.weixin.qq.com/s/cREIWIcDftXgVyTztGjV3Q)
  - [单个10x单细胞转录组的样品有多个fq文件该如何定量](https://mp.weixin.qq.com/s/UgDGwDeufaRUjEzzU90-HQ)