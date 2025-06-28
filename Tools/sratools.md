# Use sra-tools to download big data from NCBI

**Reference**: [测序数据处理 —— 数据下载](https://mp.weixin.qq.com/s/rxlfIRZ7AhdXNlDxOMys0w)
---

# [Install](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)

```shell
# Old image: stereonote-multi-gcc
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
cd ~
cp /data/work/sratoolkit.tar.gz .
tar -vxzf sratoolkit.tar.gz
export PATH=$PWD/sratoolkit.3.2.1-ubuntu64/bin:$PATH
rm sratoolkit.tar.gz
# New image: sratools
```

# Usages
