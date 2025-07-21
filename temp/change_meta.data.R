# 创建一个空的映射表
# mapping <- data.frame(cluster = character(), celltype = character(), stringsAsFactors = FALSE)

# 遍历每一行，将 cluster 列的每个值与 celltype 列的值进行匹配
# for (i in 1:nrow(df)) {
#   clusters <- df$cluster[[i]]
#   celltype <- df$celltype[i]
#   mapping <- rbind(mapping, data.frame(cluster = clusters, celltype = celltype, stringsAsFactors = FALSE))
# }
mapping <- read.csv('old_new.csv')
seu$rename0626 <- ifelse(seu$rename0626 %in% mapping$old, mapping$new[match(seu$rename0626, mapping$cluster)], seu$rename0626)
# 查看映射表
print("\n映射表：")
print(mapping)
write.csv(mapping, paste0(string_name,".csv"), row.names = FALSE) 

colnames(seu@meta.data)
seu$rename0626 <- ifelse(seu$rename0626 %in% mapping$cluster, mapping$celltype[match(seu$rename0626, mapping$cluster)], seu$rename0626)

代码解释