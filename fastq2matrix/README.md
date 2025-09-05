# fastq2matrix: Embrace shared data

---
```shell
cd /opt
# [ download file from downloads page ]
wget -O cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1757082956&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=EhdqMA~vIQizKQIIkf9QIF04dpc3bWDQZm7ehNWxui~Ne1KOh0Zbs8Vchipe7hHrCD36bSc1Cn5u5E8gkv12ILO~xcZKpLlwOmYYkkz3IJGYPXkaTcLjCVW9LXGs1iAm7cNrGoeuPd~ZV1csplqx~aXytb0Kdc-RAcNWPf~2Gr2uf18yMdEDaZ~zpQq~ZaaG09yirqM-wJ9oyi0B6aNN5xNIg9qsVb63R3vdFcy0i2bMhKH~L9~6Hw9jEyuA4E5jXO98LUCwz2JX7qtN0ioWr9wRPX4GlHqScBh4cz-tWQT3NbBclbK~XRULrSNLdHxA0LQO-xk6ADN4JfCDHiWS4w__"
tar -xzvf cellranger-9.0.1.tar.gz
export PATH=/opt/cellranger-9.0.1:$PATH
rm cellranger-9.0.1.tar.gz
# wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
# tar -vxzf sratoolkit.tar.gz
# export PATH=$PWD/sratoolkit.3.2.1-ubuntu64/bin:$PATH
conda install bioconda::sra-tools -y
export PATH=/opt/software/miniconda3/bin:$PATH
```
# Reference & Citation
- [https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)
- [https://www.10xgenomics.com/support/cn/software/cell-ranger/downloads](https://www.10xgenomics.com/support/cn/software/cell-ranger/downloads)
- [https://www.10xgenomics.com/support/cn/software/cell-ranger-arc/downloads#download-links](https://www.10xgenomics.com/support/cn/software/cell-ranger-arc/downloads#download-links)
- [单细胞转录组：下载 NCBI 数据 （SRA ）：从fastq 到10x矩阵](https://mp.weixin.qq.com/s/4tqmyMEkCtn0UJwN-l4Swg)
- [单细胞专题31| Cellranger: 10x 转录组从原始下机数据到表达矩阵](https://mp.weixin.qq.com/s/j0iZFRKaQrVV3Lq4TVzWuQ)