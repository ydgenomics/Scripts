#!/bin/bash

# 定义输入文件路径
input_file="/data/work/convert/0525test/multi_layers/H1314_dataget_Anno_rename_threelayers.hr.rds"
sctype="rds2h5ad"
#sctype="h5ad2rds"
# 获取文件后缀
file_extension="${input_file##*.}"
echo "File extension: $file_extension"


if [[ "$sctype" == "rds2h5ad" ]]; then
    if [[ "$file_extension" == "rds" ]]; then
        echo "Converting rds to h5ad..."
        /software/conda/Anaconda/bin/Rscript /script/convert_rdsAh5ad2.R \
        --input_file $input_file
    elif [[ "$file_extension" == "h5ad" ]]; then
        echo "Copying h5ad file..."
        cp "$input_file" ./
    else
        echo "Unsupported file type: $file_extension"
    fi
elif [[ "$sctype" == "h5ad2rds" ]]; then
    if [[ "$file_extension" == "h5ad" ]]; then
        echo "Converting h5ad to rds..."
        /software/conda/Anaconda/bin/Rscript /script/convert_rdsAh5ad2.R \
        --input_file $input_file
    elif [[ "$file_extension" == "rds" ]]; then
        echo "Copying rds file..."
        cp "$input_file" ./
    else
        echo "Unsupported file type: $file_extension"
    fi
else
    echo "Unsupported conversion type: $sctype"
fi
