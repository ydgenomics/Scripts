
#生成注释信息：



# 添加热图并设置背景注释色块，包含行和列注释信息
Heatmap(as.matrix(aver_dt),
    width = unit(w, "cm"), # 设置热图宽度
    height = unit(h, "cm"), # 设置热图高度
    name = "expression", # 图例名称
    col = mycol, # 热图颜色
    cluster_columns = FALSE, # 不聚类列
    cluster_rows = FALSE, # 不聚类行
    column_names_side = "top", # 列名显示在顶部
    column_names_rot = 60, # 列名旋转角度
    row_title = "Top 5 TFs", # 行标题
    rect_gp = gpar(col = "white", lwd = 1.5), # 热图边框设置
    heatmap_legend_param = list(
        legend_height = unit(2.8, "cm"), # 图例高度
        grid_width = unit(0.4, "cm"), # 图例色块宽度
        labels_gp = gpar(col = "gray20", fontsize = 10) # 图例标签样式
    ),
    top_annotation = col_anno # 添加列注释
) +
    row_anno # 添加行注释

col_anno <- HeatmapAnnotation(
  df = aver_dt,
  show_annotation_name = FALSE,
  gp = gpar(col = "white", lwd = 2),
  col = list(
    cell = c(
    "g0" = "#E41A1C", "g1" = "#377EB8", "g2" = "#4DAF4A", "g3" = "#984EA3", "g4" = "#FF7F00",
    "g5" = "#FFFF33", "g6" = "#A65628", "g7" = "#F781BF", "g8" = "#999999", "g9" = "#66C2A5",
    "g10" = "#FC8D62", "g11" = "#8DA0CB", "g12" = "#E78AC3", "g13" = "#A6D854", "g14" = "#FFD92F",
    "g15" = "#E5C494", "g16" = "#B3B3B3", "g17" = "#B15928", "g18" = "#1F78B4"
    )
  )
)

# 根据 cluster 分类给 gene 添加颜色注释
# 假设 result_df 有 cluster 和 gene 两列
cluster_colors <- c(
    "dividing cells" = "#E41A1C",
    "Meristem" = "#377EB8",
    "procambium" = "#4DAF4A",
    "Phloem" = "#984EA3",
    "Xylem" = "#FF7F00",
    "Xylem parenchyma" = "#FFFF33",
    "Sclerenchyma" = "#A65628",
    "Epidermis" = "#F781BF",
    "Large parenchyma" = "#999999",
    "Bundle sheath" = "#66C2A5"
)
row_cols <- setNames(cluster_colors[result_df$cluster], result_df$gene)
row_cols
head(result_df, n = 8)
# cluster	gene
# <chr>	<chr>
# 1	dividing cells	LOC-Os04g40940
# 2	dividing cells	LOC-Os01g16650
# 3	dividing cells	LOC-Os04g47580
# 4	dividing cells	LOC-Os08g40170
# 5	Meristem	LOC-Os07g03770
# 6	Meristem	LOC-Os03g51690
# 7	Meristem	LOC-Os01g19694
# 8	Meristem	LOC-Os10g26340

#生成注释信息：
row_anno <- rowAnnotation(foo = anno_text(rownames(aver_dt),
								location = 0,
								just = "left",
								gp = gpar(fill = row_cols,
										  col = "black",
										  fontface = 'italic'),
								width = max_text_width(rownames(dt))*1.25))


# 添加热图并设置背景注释色块，包含行和列注释信息
Heatmap(as.matrix(aver_dt),
    width = unit(w, "cm"), # 设置热图宽度
    height = unit(h, "cm"), # 设置热图高度
    name = "expression", # 图例名称
    #col = mycol, # 热图颜色
    cluster_columns = FALSE, # 不聚类列
    cluster_rows = FALSE, # 不聚类行
    column_names_side = "top", # 列名显示在顶部
    column_names_rot = 60, # 列名旋转角度
    row_title = "Top 5 TFs", # 行标题
    rect_gp = gpar(col = "white", lwd = 1.5), # 热图边框设置
    heatmap_legend_param = list(
        legend_height = unit(2.8, "cm"), # 图例高度
        grid_width = unit(0.4, "cm"), # 图例色块宽度
        labels_gp = gpar(col = "gray20", fontsize = 10) # 图例标签样式
    ),
    top_annotation = col_anno, # 添加列注释
    left_annotation = row_anno, # 添加行注释
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = unit(0.8, "cm"), height = unit(0.3, "cm"), gp = gpar(col = NA, fill = fill))
    }
)
# 添加热图并设置背景注释色块，包含行和列注释信息
pdf('test.pdf',width=12, height=24)
Heatmap(as.matrix(aver_dt),
    # width = unit(w, "cm"), # 设置热图宽度
    # height = unit(h, "cm"), # 设置热图高度
    name = "expression", # 图例名称
    #col = mycol, # 热图颜色
    cluster_columns = FALSE, # 不聚类列
    cluster_rows = FALSE, # 不聚类行
    column_names_side = "top", # 列名显示在顶部
    column_names_rot = 60, # 列名旋转角度
    row_title = "Top 5 TFs", # 行标题
    rect_gp = gpar(col = "white", lwd = 1.5), # 热图边框设置
    heatmap_legend_param = list(
        legend_height = unit(2.8, "cm"), # 图例高度
        grid_width = unit(0.4, "cm"), # 图例色块宽度
        labels_gp = gpar(col = "gray20", fontsize = 10) # 图例标签样式
    ),
    top_annotation = col_anno, # 添加列注释
    # right_annotation = row_anno, # 添加行注释
    # cell_fun = function(j, i, x, y, width, height, fill) {
    #     grid.rect(x = x, y = y, width = unit(0.8, "cm"), height = unit(0.3, "cm"), gp = gpar(col = NA, fill = fill))
    # }
)
dev.off()