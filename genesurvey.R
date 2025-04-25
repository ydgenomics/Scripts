# Date: 2025-04-01
# Demo: 

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(tibble)

# Load necessary library for command-line arguments
library(optparse)
option_list <- list(
  make_option(c("-r", "--rds"), type = "character", default = "/data/work/tomato/1.annotation/SixTime_SCT_cellannotation.rds", # nolint
              help = "Path to the RDS file", metavar = "character"),
  make_option(c("-c", "--genecsv"), type = "character", default = "/data/work/tomato/1.genesurvey/Phase_genes_ITAG4.1.csv", # nolint
              help = "Path to the gene CSV file", metavar = "character"),
  make_option(c("-a", "--assay"), type = "character", default = "SCT",
              help = "Assay to use", metavar = "character"),
  make_option(c("-o", "--order"), type = "character", default = "TM0,TM1,TM2,LTM,EFM,SIM", # nolint
              help = "Order of samples", metavar = "character"),
  make_option(c("-s", "--sample_var"), type = "character", default = "sample",
              help = "Sample variable", metavar = "character"),
  make_option(c("-g", "--group_var"), type = "character", default = "assign.ident", # nolint
              help = "Group variable", metavar = "character")
)

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list = option_list))
rds <- opt$rds
genecsv <- opt$genecsv
assay <- opt$assay
order <- opt$order
sample_var <- opt$sample_var
group_var <- opt$group_var

gene_column <- "Genename_ITAG4.1" # nolint

seu <- readRDS(rds)
DefaultAssay(seu) <- assay
Idents(seu) <- group_var
gene_names <- read.csv(genecsv, header = TRUE, stringsAsFactors = FALSE)[[gene_column]] # nolint
head(gene_names)
gene_names <- gene_names[10:12]
seu

order <- strsplit(order, ",")[[1]]
seu <- seu[, order(match(seu@meta.data[[sample_var]], order))] # change the oder of samples # nolint
seu@meta.data[[sample_var]] <- factor(seu@meta.data[[sample_var]], levels = order) # change the oder of samples # nolint

pdf(paste0("featur-vin_plots_",sample_var,"_",group_var,"_",assay,".pdf"),  width = 5*length(order), height = 5*length(gene_names)) # nolint
# options(repr.plot.width = 18, repr.plot.height = 36) # nolint
p1 <- FeaturePlot(seu, features = gene_names, split.by = sample_var, pt.size = 0.005, label = TRUE, label.size = 3, label.color = "black")  # nolint
print(p1)
p2 <- VlnPlot(seu, assay = assay, features = gene_names, group.by = group_var, split.by = sample_var, pt.size = 0.1, stack = T, flip = T, fill.by = "feature",split.plot = TRUE) # nolint
print(p2)
dev.off()

pdf(paste0("heatmap_", sample_var, "_", group_var, "_", assay, ".pdf"))
for (gene in gene_names) {
  gene_expr <- FetchData(seu, vars = gene, slot = "data")
  gene_expr_df <- data.frame(Sample = seu[[sample_var]],Group = seu[[group_var]],  Expression = gene_expr,stringsAsFactors = FALSE) # nolint
  colnames(gene_expr_df) <- c("Sample", "Group", "Expression")
  # Calculate the mean expression for each Sample and Group combination
  mean_expr_matrix <- gene_expr_df %>%
    group_by(Sample, Group) %>%
    summarise(Mean_Expression = mean(Expression, na.rm = TRUE)) %>%
    spread(key = Group, value = Mean_Expression)  # Use spread to reshape data [^49^] # nolint
  rownames(mean_expr_matrix) <- mean_expr_matrix$Sample
  mean_expr_matrix <- as.matrix(mean_expr_matrix)
  mean_expr_matrix <- mean_expr_matrix[, -1]
  mean_expr_matrix <- as.matrix(mean_expr_matrix)
  mean_expr_matrix <- apply(mean_expr_matrix, 2, function(x) {
    if (!is.numeric(x)) {
      as.numeric(as.character(x))
    } else {
      x
    }
  })
  rownames(mean_expr_matrix) <- order
  pheatmap(mean_expr_matrix,
           #scale = "column",  # 对列进行标准化
           cluster_rows = FALSE,  # 禁用行聚类
           cluster_cols = FALSE,  # 禁用列聚类
           show_rownames = TRUE,
           show_colnames = TRUE,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           main = paste0("Mean Expression Heatmap in ", gene))
  pheatmap(mean_expr_matrix,
           scale = "column",  # 对列进行标准化
           cluster_rows = FALSE,  # 禁用行聚类
           cluster_cols = FALSE,  # 禁用列聚类
           show_rownames = TRUE,
           show_colnames = TRUE,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           main = paste0("Mean Expression Heatmap in ", gene))
}
dev.off()