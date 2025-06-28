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


# What
### NCBI SRA、SRX、SRR 的含义

#### 1. SRA（Sequence Read Archive）
SRA 是 NCBI（美国国家生物技术信息中心）维护的一个公共数据库，专门用于存储高通量测序产生的原始数据，如 DNA 测序（基因组）、RNA 测序（转录组）、元基因组、表观遗传组等。SRA 数据库是国际核苷酸序列数据库合作组织（INSDC）的一部分，由 NCBI、欧洲生物信息学研究所 (EBI) 和日本 DNA 数据库 (DDBJ) 共同维护。

#### 2. SRX（Experiment）
SRX 表示一个实验（Experiment），是 SRA 数据库中的一个层级，记录了实验设计、实验平台和结果处理等信息。一个实验可以包含多个样本（Sample）和多个运行（Run）。实验是 SRA 数据库的最基本单元，一个实验信息可以同时包含多个结果集（run）。

#### 3. SRR（Run）
SRR 表示一个运行（Run），是 SRA 数据库中的一个层级，记录了测序仪运行所产生的 reads。一个 Run 包括测序序列及质量数据。Runs 是实际存储 DNA 序列数据的地方，是你从 NCBI 下载数据时需要关注的部分。

### 数据库的层级结构
SRA 数据库的数据结构分为四个层级，从高到低依次为：
- **Study（项目）**：以 SRP 开头，表示一个研究项目，可以包含多个实验（Experiment）。
- **Experiment（实验）**：以 SRX 开头，表示一个实验，包含实验设计、样本信息、测序平台等信息。
- **Sample（样本）**：以 SRS 开头，表示样本信息，包括物种信息、菌株信息、组织类型等。
- **Run（运行）**：以 SRR 开头，表示测序仪运行所产生的 reads，是实际的测序数据。
