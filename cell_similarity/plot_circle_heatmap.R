# [R绘图：环形热图绘制](https://mp.weixin.qq.com/s/bQuZraaz65Rv2CSlHi14-w)
# [ggcor 的环形热图](https://mp.weixin.qq.com/s/tVxalBWsxLn58RJkpb-PaQ)

p8 <- quickcor(celltype_NV, circular = TRUE, cluster = TRUE, grid.colour = 'white',
			   open = 45, # 缺口大小
			   # 内圈外圈比例
			   outer = 0.2, inner = 2) +
  # 单元格边框线颜色
  geom_colour(colour = 'white') +
  # 自定义填充颜色
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') +
  # 更改图例名称
  guides(fill = guide_colorbar(title = 'AUROC')) +
  # 添加聚类树
  anno_col_tree(bcols = rand_color(5)) +
  anno_row_tree(bcols = rand_color(10)) +
  # 基因名
  set_p_yaxis(bcols = rand_color(10)) +
  # 样本名
  set_p_xaxis(bcols = rand_color(5))

quickcor(celltype_NV, circular = TRUE, cluster = TRUE, grid.colour = 'white',
         open = 45, # 缺口大小
         # 内圈外圈比例
         outer = 0.5, inner = 1) +
  # 单元格边框线颜色
  geom_colour(colour = 'white') +
  # 自定义填充颜色
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') +
  # 更改图例名称
  guides(fill = guide_colorbar(title = 'log2FC')) +
  # 列注释
  anno_hc_bar(k = 2, fill = rand_color(2), pos = 'top', height = 0.3) +
  anno_hc_bar(k = 3, fill = rand_color(3), pos = 'top', height = 0.5) +
  anno_hc_bar(k = 5, fill = rand_color(5), pos = 'top', height = 0.2) +
  # 添加聚类树
  anno_col_tree(bcols = rand_color(5), height = 0.15) +
  anno_hc_bar(k = 15, fill = rand_color(15), pos = 'left', width = 0.5) +
  anno_hc_bar(k = 10, fill = rand_color(10), pos = 'left', width = 0.5) +
  anno_row_tree(bcols = rand_color(8)) +

  # 样本名
  set_p_xaxis(bcols = rand_color(5)) +
  # 基因名
  set_p_yaxis()

p7 <- quickcor(da, circular = TRUE, cluster = TRUE,
         open = 45, # 缺口大小
         # 内圈外圈比例
         outer = 0.5, inner = 1.5) +
  # 单元格边框线颜色
  geom_colour(colour = 'black') +
  # 自定义填充颜色
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') +
  # 更改图例名称
  guides(fill = guide_colorbar(title = 'log2FC')) +
  # 添加聚类树
  anno_col_tree() +
  anno_row_tree() +
  # 基因名
  set_p_yaxis() +
  # 样本名
  set_p_xaxis()