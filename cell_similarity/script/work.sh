rds="~{rds}"
output_name="~{output_name}"
batch_key="~{batch_key}"
cluster_key="~{cluster_key}"
threshold_value=0.95

/opt/conda/bin/Rscript /script/cell_similarity/jaccard_hclust.R \
--input_file $rds --output_name $output_name --batch_key $batch_key --cluster_key $cluster_key
rm Rplots.pdf

/opt/conda/bin/Rscript /script/cell_similarity/metaNeighbor.R \
--input_file $rds --output_name $output_name --batch_key $batch_key --cluster_key $cluster_key --threshold_value $threshold_value

path=$(find "$(pwd)" -maxdepth 1 -name '*_metaNeighbor.csv' -exec readlink -f {} \;)
path=$(echo "$path" | head -n 1)
echo "Path to celltype_NV of metaNeighbor output: $path"
seq='seq.txt'
slimit=0.95
python /script/cell_similarity/sanky_plot.py \
--path $path --seq $seq --slimit $slimit