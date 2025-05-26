import scanpy as sc
import os

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