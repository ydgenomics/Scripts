# date: 20250609
# image: go-figure--01
# parameters

# tsv_path="/data/work/go-figure/standard_example_input.tsv"
# max_label=15
# # run go-figure
# mkdir data
# cp /script/go-figure/ic.tsv ./data
# cp /script/go-figure/relations_full.tsv ./data
# cp /script/go-figure/go.obo ./data
# mkdir go-figure-output
# /software/miniconda/envs/go-figure/bin/python /script/go-figure/gofigure.py \
# -i $tsv_path -j standard -m $max_label -o go-figure-output
# 定义文件路径
result_name_file="../result_name.txt"
output_file="../output_standard_gofigure_input.txt"

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
  max_label=15

  # run go-figure
  mkdir "./$name/data"
  cp /data/work/go-figure/data/ic.tsv "$name/data"
  cp /data/work/go-figure/data/relations_full.tsv "$name/data"
  cp /data/work/go-figure/data/go.obo "$name/data"
  mkdir "$name/go-figure-output"
  /software/miniconda/envs/go-figure/bin/python /data/work/go-figure/gofigure.py \
    -i "$tsv_path" -j standard -m "$max_label" -o "$name/go-figure-output"
  rm -r "$name/data"

  # 自增索引
  i=$((i + 1))
done