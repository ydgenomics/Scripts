# Tilte:rice_annotation.R
# Date: 2025-04-30
#"
#cd /data/work/visual_rice_annotation/2.D3
#/software/miniconda/envs/Seurat/bin/Rscript /data/work/visual_rice_annotation/rice_annotation.R \
#--testrds_path /data/users/yangdong/yangdong_f6fd22e6e3a247ceaae97934225564ba/online/result/scPlant/test/NipLSD3_obj_after_choir.rds \
#--celljson_path /data/users/yangdong/yangdong_f6fd22e6e3a247ceaae97934225564ba/online/result/scPlant/CRA004082_Rice_celltype.vocab.meta.json \
#--data_path /data/users/yangdong/yangdong_f6fd22e6e3a247ceaae97934225564ba/online/output/rice_scplantllm_D3/combined_data.csv \
#--group_var "cell_type,singleR,sctype_classification,predict_cell,max_cell" \
#--reference_rdata "/data/users/yangdong/yangdong_f6fd22e6e3a247ceaae97934225564ba/online/result/scPlant/Reference_CRA004082_leaf_filtered_cg_pp.Rdata" \
#--db_ "/data/users/yangdong/yangdong_f6fd22e6e3a247ceaae97934225564ba/online/script/sctype/rice_leaf_marker_raw.xlsx" \
#--tissue "leaf"
#"

library(Seurat)
library(jsonlite)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(scales)
library(RColorBrewer)
library(grDevices)
library(optparse)
library(gridExtra)
library(grid)
library(gridExtra)
lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(SingleR)

option_list <- list(
    make_option(c("-t", "--testrds_path"), type = "character", default = "/data/users/yangdong/yangdong_f6fd22e6e3a247ceaae97934225564ba/online/result/scPlant/test/NipLSD1_obj_after_choir.rds", 
                            help = "Path to the testrds file", metavar = "character"),
    make_option(c("-c", "--celljson_path"), type = "character", default = "/data/users/yangdong/yangdong_f6fd22e6e3a247ceaae97934225564ba/online/result/scPlant/CRA004082_Rice_celltype.vocab.meta.json", 
                            help = "Path to the celltype vocab JSON file", metavar = "character"),
    make_option(c("-d", "--data_path"), type = "character", default = "/data/users/yangdong/yangdong_f6fd22e6e3a247ceaae97934225564ba/online/output/rice_scplantllm_D1/combined_data.csv", 
                            help = "Path to the combined data CSV file", metavar = "character"),
    make_option(c("-g", "--group_var"), type = "character", default = "cell_type,singleR,sctype_classification,predict_cell,max_cell", # nolint
                            help = "Variable name for groups"),
    make_option(c("-r", "--reference_rdata"), type = "character", default = "/data/users/yangdong/yangdong_f6fd22e6e3a247ceaae97934225564ba/online/result/scPlant/Reference_CRA004082_leaf_filtered_cg_pp.Rdata", 
                            help = "Path to the reference RData file", metavar = "character"),
    make_option(c("-s", "--db_"), type = "character", default = "/data/users/yangdong/yangdong_f6fd22e6e3a247ceaae97934225564ba/online/script/sctype/rice_leaf_marker_raw.xlsx",
                            help = "Path to the ScTypeDB file", metavar = "character"),
    make_option(c("-o", "--tissue"), type = "character", default = "leaf",
                            help = "Tissue type for ScTypeDB", metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
testrds_path <- opt$testrds_path
celljson_path <- opt$celljson_path
data_path <- opt$data_path
group_var <- opt$group_var
reference_rdata <- opt$reference_rdata
db_ <- opt$db_
tissue <- opt$tissue

group_var_list <- unlist(strsplit(group_var, ","))

#testrds_path="/data/work/scPlantLLM/data/testrds.txt"
#data <- read.csv(testrds_path, header = FALSE, sep = ",", stringsAsFactors = FALSE, col.names = c("path"))
#file_paths <- as.list(data$path)
#print(file_paths)

#if (length(file_paths) > 1) {
#    merged_data <- readRDS(file_paths[[1]])
#    for (i in 2:length(file_paths)) {
#        temp_data <- readRDS(file_paths[[i]])
#        merged_data <- merge(merged_data, temp_data)
#    }
#    print(merged_data)
#} else {
#    merged_data <- readRDS(file_paths[[1]])
#    print(merged_data)
#}
merged_data <- readRDS(testrds_path)
Idents(merged_data) <- "CHOIR_clusters_0.05"

celltype_vocab <- fromJSON(celljson_path)
data <- read.csv(data_path, header = TRUE, sep = ",", stringsAsFactors = FALSE)
vocab_df <- data.frame(
  Cell_Type_Prediction = unlist(celltype_vocab),
  Cell_Type = names(celltype_vocab),
  stringsAsFactors = FALSE
)
print(vocab_df)
data <- left_join(data, vocab_df, by = "Cell_Type_Prediction")
print(head(data))

######################################### predict_cell #########################################
# Match cell_name of merged_data 
cell_names <- colnames(merged_data)
predict_cell <- rep(NA, length(cell_names))
match_indices <- match(cell_names, data$Cell_Name)
predict_cell[!is.na(match_indices)] <- data$Cell_Type[match_indices[!is.na(match_indices)]]
merged_data$predict_cell <- predict_cell
print(head(merged_data)) # predict_cell
unique(merged_data$predict_cell)
# Visual
#merged_data <- NormalizeData(merged_data)
#merged_data <- FindVariableFeatures(merged_data, nfeatures = 3000)
#merged_data <- ScaleData(merged_data)
#merged_data <- RunPCA(merged_data, features = VariableFeatures(object = merged_data), verbose = FALSE)
#merged_data <- RunUMAP(merged_data, dims = 1:20, verbose = FALSE)

# Max_cell type prediction
cluster_counts <- table(merged_data$CHOIR_clusters_0.05)
cell_type_counts <- merged_data@meta.data %>%
  as.data.frame() %>%
  group_by(CHOIR_clusters_0.05, predict_cell) %>%
  summarise(n = n()) %>%
  group_by(CHOIR_clusters_0.05) %>%
  mutate(percent = n / sum(n) * 100)

p <- ggplot(cell_type_counts, aes(x = CHOIR_clusters_0.05, y = n, fill = predict_cell)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.3) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "Clusters",
    y = "Number of Cells",
    fill = "Predicted Cell Type",
    title = "Cell Type Distribution"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(fill = guide_legend(title = "Predicted Cell Type"))

############################################ Max_cell type prediction ############################################
merged_data$max_cell <- as.character(merged_data$CHOIR_clusters_0.05) # max_cell
merged_data@meta.data$max_cell <- as.character(merged_data@meta.data$max_cell)
unique_clusters <- unique(merged_data$CHOIR_clusters_0.05)
for (cluster in unique_clusters) {
  cluster <- as.character(cluster)  # 确保cluster是字符类型
  cell_types <- cell_type_counts[cell_type_counts$CHOIR_clusters_0.05 == cluster, ]
  max_cell_type <- as.character(cell_types[which.max(cell_types$percent), "predict_cell"])
  merged_data@meta.data$max_cell[merged_data@meta.data$max_cell == cluster] <- max_cell_type
  merged_data$max_cell <- merged_data@meta.data$max_cell
  print(paste("Cluster:", cluster, "=> Max Cell Type:", max_cell_type))
}
unique(merged_data$max_cell)

results <- lapply(unique(merged_data$CHOIR_clusters_0.05), function(cluster) {
  cluster_cells <- cluster_counts[as.character(cluster)]
  cell_types <- cell_type_counts[cell_type_counts$CHOIR_clusters_0.05 == cluster, ]
  cell_type_str <- paste0(
    cell_types$predict_cell,
    " ", round(cell_types$percent, 2), "%",
    collapse = ", "
  )
  max_cell_type <- cell_types[which.max(cell_types$percent), "predict_cell"]
  paste0(
    "cluster ", cluster, ": ", cluster_cells, " cell; ","probability ", max_cell_type, "; ", cell_type_str
  )
})
writeLines(unlist(results), "scplantllm_predict_summary.txt")
#saveRDS(merged_data, "merged_data.rds")

##################################### SingleR #################################
load(reference_rdata)
ref_sce
testdata <- GetAssayData(merged_data, slot = "data")
common_genes <- intersect(rownames(testdata), rownames(ref_sce)) # 找到两个数据集中都存在的基因名
num_common_genes <- length(common_genes) # 计算重叠的基因数目
print(paste0("The common gene number of Reference and Test data: ", num_common_genes)) # 输出重叠的基因数目 #"The common gene number of Reference and Test data: 26025"
pred <- SingleR(test = testdata, ref = ref_sce, labels = ref_sce$Type)
as.data.frame(table(pred$labels))
merged_data$singleR=pred$labels
unique(merged_data$singleR)
#p <- DimPlot(merged_data, reduction = "umap", repel = TRUE, group.by = "singleR", label = TRUE)

################################# sctype_classfication #################################
### Cell type assignment
# load gene set preparation function
gene_sets_prepare <- function(path_to_db_file, cell_type){
  cell_markers = openxlsx::read.xlsx(path_to_db_file)
  cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
  list(gs_positive = gs, gs_negative = gs2)
}

# load cell type annotation function
sctype_score <- function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...){
  
  # check input matrix
  if(!is.matrix(scRNAseqData)){
    warning("scRNAseqData doesn't seem to be a matrix")
  } else {
    if(sum(dim(scRNAseqData))==0){
       warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
    }
  }
  
  # marker sensitivity
  marker_stat = sort(table(unlist(gs)), decreasing = T); 
  marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
                                      gene_ = names(marker_stat), stringsAsFactors = !1)

  # convert gene names to Uppercase
  #if(gene_names_to_uppercase){
  #  rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
  #}
  
  # subselect genes only found in data
  names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
  gs = lapply(1:length(gs), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  gs2 = lapply(1:length(gs2), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]
  
  # z-scale if not
  if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData
  
  # multiple by marker sensitivity
  for(jj in 1:nrow(cell_markers_genes_score)){
    Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }
  
  # subselect only with marker genes
  Z = Z[unique(c(unlist(gs),unlist(gs2))), ]
  
  # combine scores
  es = do.call("rbind", lapply(names(gs), function(gss_){ 
    sapply(1:ncol(Z), function(j) {
      gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
      sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
      if(is.na(sum_t2)){
        sum_t2 = 0;
      }
      sum_t1 + sum_t2
    })
  })) 
  
  dimnames(es) = list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows
 
  es.max
}

# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)
str(gs_list)

pbmc <- merged_data

# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(pbmc[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# extract scaled scRNA-seq matrix
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(pbmc[["RNA"]]$scale.data) else as.matrix(pbmc[["RNA"]]@scale.data)

# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
str(es.max)

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. For raw (unscaled) count matrix set scaled = FALSE
# When using Seurat, we use "RNA" slot with 'scale.data' by default. Please change "RNA" to "SCT" for sctransform-normalized data,
# or to "integrated" for joint dataset analysis. To apply sctype with unscaled data, use e.g. pbmc[["RNA"]]$counts or pbmc[["RNA"]]@counts, with scaled set to FALSE.

# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(pbmc@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(pbmc@meta.data[pbmc@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(pbmc@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3])

pbmc@meta.data$sctype_classification = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  pbmc@meta.data$sctype_classification[pbmc@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

#p3 <- DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification')        


#source("https://raw.githubusercontent.com/kris-nader/sc-type/master/R/sctype_wrapper.R"); 
#pbmc <- run_sctype(pbmc,known_tissue_type="Immune system",custom_marker_file="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",name="sctype_classification",plot=TRUE)


# prepare edges
cL_resutls <- cL_resutls[order(cL_resutls$cluster),]; edges = cL_resutls; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL

# prepare nodes
nodes_lvl1 <- sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
#ccolss <- c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")

# 获取 Idents 的唯一值数量
num_idents <- length(unique(Idents(merged_data)))
# 预定义的颜色方案
predefined_colors <- c("#5f75ae", "#92bbb8", "#64a841", "#e5486e", "#de8e06", "#eccf5a", 
                       "#b5aa0f", "#e4b680", "#7ba39d", "#b15928", "#ffff99", 
                       "#6a3d9a", "#cab2d6", "#ff7f00", "#fdbf6f", "#e31a1c", 
                       "#fb9a99", "#33a02c", "#b2df8a", "#1f78b4", "#a6cee3")
# 如果唯一值数量大于 21，生成随机颜色
if (num_idents > 21) {
  # 获取所有可用的颜色
  all_colors <- colors()
  # 从可用颜色中随机选择不重复的颜色
  random_colors <- sample(all_colors, num_idents, replace = FALSE)
  # 将随机颜色存储到 ccolss
  ccolss <- random_colors
} else {
  # 如果唯一值数量小于或等于 21，使用预定义的颜色方案
  ccolss <- predefined_colors[1:num_idents]
}
# 查看生成的颜色方案
print(ccolss)
               
for (i in 1:length(unique(cL_resutls$cluster))){
  dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
}
nodes <- rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
files_db <- openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

mygraph <- graph_from_data_frame(edges, vertices=nodes)

# Make the graph
gggr <- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")
  
p4 <- DimPlot(pbmc, reduction = "CHOIR_P0_reduction_UMAP", label = TRUE, repel = TRUE, cols = ccolss)+ gggr


merged_data <- pbmc

saveRDS(merged_data, "merged_annotated.rds")
########################### Plotting ###########################
create_percent_plot <- function(objs, mycolor_string, count_threshold = 250, sample_var = "sample", group_var = "assign.ident", percentage_threshold = 1) { # nolint
  mycolor <- unlist(strsplit(mycolor_string, ","))
  num_groups <- length(unique(objs@meta.data[[group_var]]))
  if (num_groups > length(mycolor)) {
      mycolor <- rainbow(num_groups)
  }
  fill_scale <- scale_fill_manual(values = mycolor)
  meta.data <- select(objs@meta.data, !!sym(sample_var),!!sym(group_var)) %>%
    group_by(!!sym(sample_var), !!sym(group_var)) %>%
    summarise(count = n()) %>%
    arrange(!!sym(sample_var), desc(!!sym(group_var))) %>%
    mutate(cumsum = cumsum(count), total_count = sum(count), percentage = count / total_count * 100, adjusted_percentage = percentage * 10) # nolint
  p1 <- ggplot(data = meta.data,aes(x=!!sym(sample_var),y=count,fill=!!sym(group_var)))+ # nolint
    geom_col()+
    geom_text(
      data = filter(meta.data, count > 250),
      aes(x=sample, y=cumsum - 0.5*count, label=count),color = "white") +
    fill_scale +
    theme_cowplot() +
    theme(axis.title.x = element_blank())
  # p2
  p2 <- ggplot(data = meta.data, aes(x = !!sym(sample_var), y = percentage, fill = !!sym(group_var))) + # nolint
      geom_col(width = 0.5) +
      geom_text(
              data = filter(meta.data, percentage > 1),
              aes(x = sample, y = percentage, label = paste0(round(percentage, 1), "%")), # nolint
              color = "white",
              position = position_stack(vjust = 0.5)
      ) +
      fill_scale +
      theme_cowplot() +
      theme(axis.title.x = element_blank())
  p <- p1 + p2 + plot_layout(ncol = 2)
  # p3
  p3 <- ggplot(meta.data, aes(x = !!sym(sample_var), y = percentage, fill = !!sym(group_var))) + # nolint
    geom_bar(stat = "identity", position = "fill") +
    geom_text(aes(label = paste0(round(percentage, 1), ""), group = !!sym(group_var)), # nolint
              position = position_fill(vjust = 0.5), size = 2.5, color = "black") + # nolint
    coord_flip() +
    theme_bw() +
    ylab("") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.title = element_text(size = 7.82, face = "bold"),
          axis.text = element_text(size = 7.82, color = 'black'),
          legend.text = element_text(size = 7.82),
          plot.title = element_text(size = 7.82, face = "bold"),
          axis.line = element_line(color = 'black'),
          legend.title = element_text(size = 7.82)) +
    scale_y_continuous(position = "right", expand = c(0, 0)) +
    fill_scale
  return(list(p = p, p3 = p3))
}

create_roe_plot <- function(objs, group_var="assign.ident", sample_var="sample", mycolor_string) { # nolint
  mycolor <- unlist(strsplit(mycolor_string, ","))
  num_groups <- length(unique(objs@meta.data[[group_var]]))
  if (num_groups > length(mycolor)) {
      mycolor <- rainbow(num_groups)
  }
  Idents(objs) <- group_var
  tbl <- table(objs@meta.data[[group_var]], objs@meta.data[[sample_var]])
  res <- chisq.test(tbl)
  roe <- tbl / res$expected
  roe_df <- as.data.frame(roe)
  roe_df <- roe_df %>%
    mutate(Cell = Var1, sample = Var2, roe = Freq) %>%
    select(-Var1, -Var2, -Freq)
  pic <- roe_df %>%
    ggplot(aes(Cell, roe, fill = sample)) +
    geom_bar(stat = "identity", position = "dodge") +
    ylab("Ro/E") +
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    geom_hline(aes(yintercept = 1), colour = "#990000", linetype = "dashed") +
    scale_fill_manual(values = mycolor)
  return(pic)
}

#color_string <- paste(ccolss, collapse = ",")
mycolor_string <- "#1f77b4,#ff7f0e,#279e68,#d62728,#aa40fc,#8c564b,#e377c2,#b5bd61,#17becf,#aec7e8,#ffbb78"
#mycolor_string <- color_string
pdf("Cell_count.pdf", width = 6*length(unique(merged_data$sample)), height = 8)
for (group_var in group_var_list) {
  percent_plots <- create_percent_plot(merged_data, group_var = group_var, mycolor_string)
  roe_plot <- create_roe_plot(merged_data, group_var = group_var, sample_var = "sample", mycolor_string = "lightblue,purple,orange,green,red")
  if (group_var == "sctype_classification"){
    print(p4)
  }
  if (is.list(percent_plots)) {
    for (plot_name in names(percent_plots)) {
      print(percent_plots[[plot_name]])
    }
  }  
  print(roe_plot)
}
print(p)
dev.off()

pdf("Cell_dim.pdf", width = 10, height = 8)
for (group_var in group_var_list) {
  p0 <- DimPlot(merged_data, reduction = "CHOIR_P0_reduction_UMAP", label = TRUE, pt.size = 0.5, group.by = group_var)
  print(p0)
  if (group_var == "sctype_classification"){
    print(p4)
  }
}
#p1 <- DimPlot(merged_data, reduction = "CHOIR_P0_reduction_UMAP", label = TRUE, pt.size = 0.5, group.by = "CHOIR_clusters_0.05")
dev.off()