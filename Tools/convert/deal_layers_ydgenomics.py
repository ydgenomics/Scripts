import scanpy as sc
import os
import re
import argparse

parser = argparse.ArgumentParser(description="Process h5ad layers.")
parser.add_argument('--input_path', type=str, default="/data/work/0.peanut/annotation/three_layers/H1314_dataget_Anno_rename_threelayers.h5ad", help='Input h5ad file path')
parser.add_argument('--sctype', type=str, default="h5ad", help='Comma-separated list of assays')
args = parser.parse_args()
input_path = args.input_path
sctype = args.sctype

def get_single_layer_h5ad(input_path):
    # 读取原始 .h5ad 文件
    adata = sc.read_h5ad(input_path)
    # 检查 adata 是否存在 .layers 属性
    if not hasattr(adata, 'layers') or 'counts' not in adata.layers:
        # 如果 .layers 属性不存在或 'counts' 层不存在，则创建 'counts' 层
        adata.layers['counts'] = adata.X.copy()
    print(adata.layers.keys())
    # 提取文件名（不包括扩展名）
    basename = os.path.basename(input_path)
    basename_without_ext = os.path.splitext(basename)[0]
    # 保存 'counts' 层的文件路径
    counts_path = os.path.abspath(f"{basename_without_ext}_counts.h5ad")
    adata.X = adata.layers['counts'].copy()
    adata.write(counts_path, compression="gzip")
    print(f"Layer: counts, Saved to: {counts_path}")
    # 获取所有层的名称
    layer_names = adata.layers.keys()
    # 存储所有保存的文件路径
    saved_paths = [counts_path]
    saved_layers = ['counts']
    # 遍历所有层
    for layer_name in layer_names:
        if layer_name != 'counts':  # 已保存 'counts' 层，跳过
            # 创建一个新的 AnnData 对象
            new_adata = sc.AnnData(X=adata.layers[layer_name].copy(), obs=adata.obs, var=adata.var)
            output_filename = f"{basename_without_ext}_{layer_name}.h5ad"
            output_path = os.path.abspath(output_filename)
            new_adata.write(output_path, compression="gzip")
            print(f"Layer: {layer_name}, Saved to: {output_path}")
            saved_paths.append(output_path)
            saved_layers.append(layer_name)
    # 将 saved_paths 保存为以逗号分隔的文本文件
    with open('saved_paths.txt', 'w') as f:
        f.write(','.join(saved_paths))

    # 将 saved_layers 保存为以逗号分隔的文本文件
    with open('saved_layers.txt', 'w') as f:
        f.write(','.join(saved_layers))

    print("Saved paths and layers to text files.")
    return saved_layers, saved_paths

def get_multi_layers_h5ad(input_path, assays, output_path):
    assays = ["counts" if assay == "RNA" else assay for assay in assays]
    adata = sc.read_h5ad(input_path[0])
    adata.layers[assays[0]] = adata.X.copy()
    for i in range(1, len(input_path)):
        adata2 = sc.read_h5ad(input_path[i])
        adata.layers[assays[i]] = adata2.X.copy()
    print(adata)
    adata.X = adata.layers['counts'].copy()
    adata.write(output_path, compression="gzip")
    for file_path in input_path:
        if os.path.exists(file_path):  # 检查文件是否存在
            os.remove(file_path)       # 删除文件
            print(f"Deleted file: {file_path}")
        else:
            print(f"File does not exist, could not delete: {file_path}")

if sctype == "h5ad":
    get_single_layer_h5ad(input_path)
else:
    file_name = os.path.basename(input_path)
    output_path = re.sub(r'\.rds$', '.rh.h5ad', file_name)
    print("output_path:", output_path)
    with open('saved_layers.txt', 'r') as file:
        saved_layers_str = file.read().strip()  # 读取文件内容并去除首尾空格
    saved_layers = saved_layers_str.split(','); print(saved_layers)
    with open('saved_paths.txt', 'r') as file:
        saved_paths_str = file.read().strip()  # 读取文件内容并去除首尾空格
    saved_paths = saved_paths_str.split(','); print(saved_paths)
    get_multi_layers_h5ad(saved_paths, saved_layers, output_path)