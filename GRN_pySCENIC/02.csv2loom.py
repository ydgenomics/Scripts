# Title: 02.csv2loom.py
# Author: ydgenomics
# Date: 20250520
# Image: 
# Description: This script is the preparation of runing pySCENIC.
# Pre-requiry:
# params: rds_path

import os,sys
import loompy as lp
import numpy as np
import scanpy as sc
import argparse

parser = argparse.ArgumentParser(description="Example script")
parser.add_argument("-n", "--name", type=str, help="Your name")
parser.add_argument("-a", "--age", type=int, help="Your age")

scenic_csv_path="scenic.data.csv"
scenic_loom_path="scenic.loom"

x=sc.read_csv(fiscenic_csv_pathle1);
row_attrs={"Gene":np.array(x.var_names),};
col_attrs={"CellID":np.array(x.obs_names)};
lp.create(scenic_loom_path,x.X.transpose(),row_attrs,col_attrs);