# Title: 02.csv2loom.py
# Author: ydgenomics
# Date: 2025-04-22
# Description: This script is the preparation of runing pySCENIC.
# Pre-requiry:
# params: rds_path

import os,sys
import loompy as lp
import numpy as np
import scanpy as sc

scenic_csv_path="scenic.data.csv"
scenic_loom_path="scenic.loom"

x=sc.read_csv(fiscenic_csv_pathle1);
row_attrs={"Gene":np.array(x.var_names),};
col_attrs={"CellID":np.array(x.obs_names)};
lp.create(scenic_loom_path,x.X.transpose(),row_attrs,col_attrs);

# step1 grn
!pyscenic grn --num_workers 30 --output grn.tsv --method grnboost2 scenic.loom /data/users/lili10/lili10_642e569efa3b4d56a57481c396194c66/online/input/pySCENIC/TF_gene_PlantTFDB.txt