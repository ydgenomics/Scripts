# date: 20250605 # image: go-figure
# input: result_txt.txt, result_name.txt
/opt/software/R/bin/Rscript /data/work/go-figure/deal_enrich_txt.R result_txt.txt

mkdir gofigure_result
cd gofigure_result

result_name_file="../result_name.txt"
output_file="../output_standard_gofigure_input.txt"
max_label=15

# 读取文件内容并按逗号分割
IFS=',' read -r -a names <<< "$(cat "$result_name_file")"
IFS=',' read -r -a outputs <<< "$(cat "$output_file")"

# 获取数组长度
len_names=${#names[@]}
len_outputs=${#outputs[@]}

# 确保两个数组长度一致
if [ "$len_names" -ne "$len_outputs" ]; then
  echo "Error: The number of names and outputs does not match."
  exit 1
fi

# 循环处理每个名称和对应的输出文件
i=0
while [ $i -lt $len_names ]; do
  name=${names[$i]}
  output=${outputs[$i]}

  echo "Processing $name with output file $output"
  mkdir "$name"
  tsv_path="$output"

  # run go-figure
  mkdir "$name"
  cp /data/work/go-figure/data/ic.tsv "data"
  cp /data/work/go-figure/data/relations_full.tsv "data"
  cp /data/work/go-figure/data/go.obo "data"
  /software/miniconda/envs/go-figure/bin/python /data/work/go-figure/gofigure.py \
    -i "$tsv_path" -j standard -m "$max_label" -o "$name"

  i=$((i + 1))
done
rm -r "data"