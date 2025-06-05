# date: 20250605
# image: go-figure--01
# parameters
tsv_path="/data/work/go-figure/standard_example_input.tsv"
max_label=15
# run go-figure
mkdir data
cp /script/go-figure/ic.tsv ./data
cp /script/go-figure/relations_full.tsv ./data
cp /script/go-figure/go.obo ./data
mkdir go-figure-output
/software/miniconda/envs/go-figure/bin/python /script/go-figure/gofigure.py \
-i $tsv_path -j standard -m $max_label -o go-figure-output