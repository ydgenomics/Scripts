# fastq2matrix: Embrace shared data

---
- **Image**
  - fastq2matrix

```shell
cd /opt
# [ download file from downloads page ]
wget -O cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1757082956&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=EhdqMA~vIQizKQIIkf9QIF04dpc3bWDQZm7ehNWxui~Ne1KOh0Zbs8Vchipe7hHrCD36bSc1Cn5u5E8gkv12ILO~xcZKpLlwOmYYkkz3IJGYPXkaTcLjCVW9LXGs1iAm7cNrGoeuPd~ZV1csplqx~aXytb0Kdc-RAcNWPf~2Gr2uf18yMdEDaZ~zpQq~ZaaG09yirqM-wJ9oyi0B6aNN5xNIg9qsVb63R3vdFcy0i2bMhKH~L9~6Hw9jEyuA4E5jXO98LUCwz2JX7qtN0ioWr9wRPX4GlHqScBh4cz-tWQT3NbBclbK~XRULrSNLdHxA0LQO-xk6ADN4JfCDHiWS4w__"
tar -xzvf cellranger-9.0.1.tar.gz
export PATH=/opt/cellranger-9.0.1:$PATH # /opt/cellranger-9.0.1/bin/cellranger
rm cellranger-9.0.1.tar.gz
# wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
# tar -vxzf sratoolkit.tar.gz
# export PATH=$PWD/sratoolkit.3.2.1-ubuntu64/bin:$PATH
conda install bioconda::sra-tools -y
export PATH=/opt/software/miniconda3/bin:$PATH # 
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
prefetch SRR31330969 --max-size 120G

# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/990/345/GCF_007990345.1_Gossypium_hirsutum_v2.1/GCF_007990345.1_Gossypium_hirsutum_v2.1_genomic.gtf.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/990/345/GCF_007990345.1_Gossypium_hirsutum_v2.1/GCF_007990345.1_Gossypium_hirsutum_v2.1_protein.faa.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/990/345/GCF_007990345.1_Gossypium_hirsutum_v2.1/GCF_007990345.1_Gossypium_hirsutum_v2.1_genomic.fna.gz
```

复现过程找到原文的参考基因组很困难，通过NCBI的refgenome来找，你们的基因名都是LOC命名方式，这些只是标识符，从标识符到相比更具意义的基因名有不少困难

  
</details>

---
# Reference & Citation
- [https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)
- [https://www.10xgenomics.com/support/cn/software/cell-ranger/downloads](https://www.10xgenomics.com/support/cn/software/cell-ranger/downloads)
- [https://www.10xgenomics.com/support/cn/software/cell-ranger-arc/downloads#download-links](https://www.10xgenomics.com/support/cn/software/cell-ranger-arc/downloads#download-links)
- [单细胞转录组：下载 NCBI 数据 （SRA ）：从fastq 到10x矩阵](https://mp.weixin.qq.com/s/4tqmyMEkCtn0UJwN-l4Swg)
- [单细胞专题31| Cellranger: 10x 转录组从原始下机数据到表达矩阵](https://mp.weixin.qq.com/s/j0iZFRKaQrVV3Lq4TVzWuQ)