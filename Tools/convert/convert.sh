#!/bin/bash
# input_file="/data/input/Files/example/Arabidopsis.h5ad"
# layers="all"
# tools="single_convert"

input_file=$1
tools=$2
layers=$3

ext="${input_file##*.}"
echo "input file extension is: $ext"
if [ "$tools" == "single_convert" ]; then
    echo "Converting in single_convert..."
    /software/conda/Anaconda/bin/Rscript /script/convert_format2.R --input_file $input_file
    cp "$input_file" ./
elif [ "$tools" == "multi_convert" ]; then
    if [ "$ext" == "rds" ]; then
        echo "Converting rds to h5ad..."
        /software/conda/Anaconda/bin/Rscript /script/convert_rdsAh5ad2.R --input_file $input_file --layers $layers
        python /script/deal_layers_ydgenomics.py --input_path $input_file --sctype $ext
        echo "Copying rds file..."
        cp "$input_file" ./
    elif [ "$ext" == "h5ad" ]; then
        echo "Converting h5ad to rds..."
        python /script/deal_layers_ydgenomics.py --input_path $input_file --sctype $ext
        /software/conda/Anaconda/bin/Rscript /script/convert_rdsAh5ad2.R --input_file $input_file --layers $layers
        echo "Copying h5ad file..."
        cp "$input_file" ./
    else
        echo "Error: Unsupported file extension '$ext'. Only 'rds' and 'h5ad' are supported."
    fi
else
    echo "Error: Unsupported tools "$tools". Only 'single_convert' and 'multi_convert' are supported."
fi
