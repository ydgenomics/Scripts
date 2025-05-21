# Author: lili, ydgenomics
# Date: 20250521
# Prepore files: input.rds; TF_gene_PlantTFDB.txt; ITAG4.1.regions_vs_motifs.rankings.feather; ITAG4.1_MOTIF_PlantTFDB.tbl
# Output files: scenic.data.csv; scenic.loom; grn.tsv; ctx.csv; aucell.loom
# What are .feather and .tbl files?

Rscript 01.csv.R \
--input_rds /data/users/yuantingting/yuantingting_fd18e4ec5adc46dd89e2e5c607e9e7e3/online/output/PlantPhone/Seurat/Merge/EFM/RNA_T_0.5/Annotated_EFM_RNA_T_0.5.rds \
--output_csv scenic.data.csv # Image: GRN-SCENIC-R--03 /opt/conda/bin/R

python 02.csv2loom.py \
--input_csv scenic.data.csv \
--output_loom scenic.loom # Image:GRN-SCENIC-R--03 /opt/conda/bin/python

# step1 grn
pyscenic grn \
--num_workers 30 \
--output grn.tsv \
--method grnboost2 scenic.loom \
/data/users/lili10/lili10_642e569efa3b4d56a57481c396194c66/online/input/pySCENIC/TF_gene_PlantTFDB.txt # Image: pySCENIC

# step2 ctx
pyscenic ctx \
grn.tsv \
/data/users/lili10/online/input/ITAG4.1.regions_vs_motifs.rankings.feather \
--annotations_fname /data/users/lili10/lili10_642e569efa3b4d56a57481c396194c66/online/input/pySCENIC/ITAG4.1_MOTIF_PlantTFDB.tbl \
--expression_mtx_fname scenic.loom \
--mode "dask_multiprocessing" \
--output ctx.csv \
--num_workers 30 \
--mask_dropouts # Image: pySCENIC

# step3 AUCell
pyscenic aucell \
scenic.loom \
ctx.csv \
--output aucell.loom \
--num_workers 30 # Image: pySCENIC

Rscript 03.plot.R \
--aucell_loom /data/work/output/pySCENIC/EFM/aucell.loom \
--input_rds /data/users/yuantingting/yuantingting_fd18e4ec5adc46dd89e2e5c607e9e7e3/online/output/PlantPhone/Seurat/Merge/EFM/RNA_T_0.5/Annotated_EFM_RNA_T_0.5.rds \
--cluster_key assign.ident # Image: GRN-SCENIC-R--03 /opt/conda/bin/R