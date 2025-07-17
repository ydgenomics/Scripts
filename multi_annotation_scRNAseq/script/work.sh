mkdir /data/work/sctype/plot
cd /data/work/sctype/plot

input_rds="/data/work/script/sctype0614/output/day3/NipLSD3_anno_merged_data_obj_after_choir_sctype.rds"
markers_csv="/data/work/script/sctype0614/rice_leaf_marker0614.csv"
cluster_color_csv="/data/work/sctype/cluster_color.csv"
cell_type="leaf"
cluster_key="sctype_resolution0.8"
reduction_key="umap"

/opt/conda/bin/Rscript /data/work/sctype/visual_cg.R \
--input_rds $input_rds --markers_csv $markers_csv \
--cluster_color_csv $cluster_color_csv --cell_type $cell_type \
--cluster_key $cluster_key --reduction_key $reduction_key


cd /data/work/sctype/anno
input_rds="/data/work/script/sctype0614/output/day3/NipLSD3_anno_merged_data_obj_after_choir_sctype.rds"
markers_csv="/data/work/script/sctype0614/rice_leaf_marker0614.csv"
cluster_color_csv="/data/work/sctype/cluster_color.csv"
cell_type="leaf"
cluster_key="resolution0.8"
reduction_key="umap"
n_circle=5

/software/miniconda/envs/Seurat/bin/Rscript /data/work/sctype/anno_sctype.R \
--input_query_rds $input_rds --input_marker_csv $markers_csv \
--tissue $cell_type --n_circle $n_circle \
--cluster_key $cluster_key --umap_name $reduction_key