library(Seurat)
library(tidyverse)

# 加载 Seurat 对象
sc <- readRDS("/data/work/input/Cer_test_BBKNNR_integrated.rds")

# 设置抽样参数
sampling_ratio <- 0.1
sampling_cells <- 200

# 提取 Seurat 的元数据
df_meta <- sc@meta.data

# 进行分组抽样
sampled_meta <- df_meta %>%
  rownames_to_column(var = "barcode") %>%
  group_by(celltype) %>%  # 假设元数据中有 "celltype" 列
  sample_n(max(sampling_cells, round(n() * sampling_ratio))) %>%
  ungroup()

# 提取抽样的细胞名
sampled_cells <- as.data.frame(sampled_meta) %>%
  pull(barcode)

# 验证抽样细胞是否存在于 Seurat 对象中
sampled_cells <- intersect(sampled_cells, colnames(sc))

# 根据抽样细胞创建新的 Seurat 对象
sc_sampled <- subset(sc, cells = sampled_cells)
saveRDS(sc_sampled, "/data/work/input/Cer_0.1_sampled.rds")