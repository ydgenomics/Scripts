[2025(BMC bioinformatics)_Gene2role: a role-based gene embedding method for comparative analysis of signed gene regulatory networks]()

# 环境配置 
> gene2role /opt/software/miniconda3/envs/gene2role/bin/
```shell
source /opt/software/miniconda3/bin/activate
conda create -n gene2role r-base=4.3 python=3.12 -y
conda activate gene2role
conda install conda-forge::r-seurat -y

pip install futures
conda install conda-forge::gensim -y #gensim 4.3.2 would require python >=3.12,<3.13.0a0 , which can be installed;
conda install conda-forge::fastdtw -y

conda install conda-forge::pandas -y
conda install conda-forge::scipy -y
conda install conda-forge::matplotlib -y

```
```R
library(Seurat)
library(data.table)
```
```python
# conda install python=3.13
import subprocess
import argparse
import os
import pandas as pd
from logging import getLogger, INFO, DEBUG
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed # pip install futures
import multiprocessing
import copy
```