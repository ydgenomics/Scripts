#!/bin/bash

# 定义输入文件路径
input_file="/data/work/convert/0525test/multi_layers/H1314_dataget_Anno_rename_threelayers.hr.rds"
sctype="rds2h5ad"
#sctype="h5ad2rds"
# 获取文件后缀
ext="${input_file##*.}"
echo "input file extension is: $ext"

if [[ "$sctype" == "rds2h5ad" ]]; then
    if [[ "$ext" == "rds" ]]; then
        echo "Converting rds to h5ad..."
        /software/conda/Anaconda/bin/Rscript /script/convert_rdsAh5ad2.R --input_file $input_file
        python /script/deal_layers_ydgenomics.py --input_path $input_file --sctype $ext
    elif [[ "$ext" == "h5ad" ]]; then
        echo "Copying h5ad file..."
        cp "$input_file" ./
    else
        echo "Unsupported file type: $file_extension"
    fi
elif [[ "$sctype" == "h5ad2rds" ]]; then
    if [[ "$ext" == "h5ad" ]]; then
        echo "Converting h5ad to rds..."
        python /script/deal_layers_ydgenomics.py --input_path $input_file --sctype $ext
        /software/script/Rscript /data/work/convert/0528test/convert_rdsAh5ad2.R --input_file $input_file 
    elif [[ "$ext" == "rds" ]]; then
        echo "Copying rds file..."
        cp "$input_file" ./
    else
        echo "Unsupported file type: $file_extension"
    fi
else
    echo "Error: Unsupported command '$sctype'. Only 'rds2h5ad' and 'h5ad2rds' are supported."
fi
