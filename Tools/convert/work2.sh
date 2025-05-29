#!/bin/bash
input_file="/data/work/0.peanut/annotation/three_layers/H1314_dataget_Anno_rename_threelayers.h5ad"
#input_file="/data/work/multi_anno/AT_root/H1314_dataget_Anno_rename_threelayers.hr.rds"
#input_file="/data/work/multi_anno/AT_root_SRP273996.rds"
#input_file="/data/work/multi_anno/AT_root/AT_root_SRP273996.rh.h5ad"
ext="${input_file##*.}"
echo "input file extension is: $ext"
if [ "$ext" == "rds" ]; then
    echo "Converting rds to h5ad..."
    /software/conda/Anaconda/bin/Rscript /script/convert_rdsAh5ad2.R --input_file $input_file --layers all
    python /script/deal_layers_ydgenomics.py --input_path $input_file --sctype $ext
    echo "Copying rds file..."
    cp "$input_file" ./
elif [ "$ext" == "h5ad" ]; then
    echo "Converting h5ad to rds..."
    python /script/deal_layers_ydgenomics.py --input_path $input_file --sctype $ext
    /software/conda/Anaconda/bin/Rscript /script/convert_rdsAh5ad2.R --input_file $input_file --layers all
    echo "Copying h5ad file..."
    cp "$input_file" ./
else
    echo "Error: Unsupported file extension '$ext'. Only 'rds' and 'h5ad' are supported."
fi