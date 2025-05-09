# Author: ydgenomics(L-plantphone) # nolint
# Date: 2025-03-28
# Description: This script generates UMAP plots, percentage plots, and Ro/E plots from a Seurat object.  # nolint
# Usage: Rscript cellsurvey.R -o <objs_path> -s <sample_var> -g <group_var> -c <mycolor_string> # nolint
# Example: Rscript cellsurvey.R -o /data/work/plantphone_test/output/_cellannotation.rds -s sample -g assign.ident -c "#1f77b4,#ff7f0e,#279e68,#d62728,#aa40fc,#8c564b,#e377c2,#b5bd61,#17becf,#aec7e8,#ffbb78" # nolint
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(scales)
library(RColorBrewer)
library(ggsci)
library(grDevices)

# Load necessary library for command-line arguments
library(optparse)
option_list <- list(
    make_option(c("-o", "--objs"), type = "character", default = "/data/work/plantphone_test/output4/TM0_TM1_cellannotation.rds", # nolint
                            help = "Path to the Seurat object RDS file"), # nolint
    make_option(c("-s", "--sample_var"), type = "character", default = "sample", # nolint
                            help = "Variable name for sample"), # nolint
    make_option(c("-g", "--group_var"), type = "character", default = "assign.ident,seurat_clusters", # nolint
                            help = "Variable name for group1"),  # nolint
    make_option(c("-c", "--mycolor_string"), type = "character",  # nolint
                            default = "#1f77b4,#ff7f0e,#279e68,#d62728,#aa40fc,#8c564b,#e377c2,#b5bd61,#17becf,#aec7e8,#ffbb78", # nolint
                            help = "Comma-separated color string"), # nolint
    make_option(c("-n", "--name"), type = "character", default = "TM0_TM1", # nolint
                            help = "Name for the output files"), # nolint
    make_option(c("-a", "--assay"), type = "character", default = "RNA", help = "Process assay") # nolint
)
opt <- parse_args(OptionParser(option_list = option_list)) # Parse command-line arguments # nolint
objs_path <- opt$objs
sample_var <- opt$sample_var
group_var <- opt$group_var
mycolor_string <- opt$mycolor_string
name <- opt$name # nolint
assay <- opt$assay

# Load Seurat object
objs <- readRDS(objs_path)
objs <- objs[, order(match(objs@meta.data$sample, c("TM0", "TM1", "TM2", "LTM", "EFM", "SIM")))] # change the oder of samples
objs@meta.data$sample <- factor(objs@meta.data$sample, levels = c("TM0", "TM1", "TM2", "LTM", "EFM", "SIM")) # change the oder of samples

DefaultAssay(objs) <- assay

#group_var <- "assign.ident,seurat_clusters"
group_var_list <- unlist(strsplit(group_var, ","))
print(group_var_list)

create_umap_plots <- function(objs, mycolor_string, sample_var = "sample", group_var = "assign.ident", reduction = "umap", label_size = 6) { # nolint
  mycolor <- unlist(strsplit(mycolor_string, ","))
  num_groups <- length(unique(objs@meta.data[[group_var]]))
  if (num_groups > length(mycolor)) {
      mycolor <- rainbow(num_groups)
  }
  pic1 <- DimPlot(objs, group.by = group_var, label = TRUE, pt.size = 2.0, reduction = reduction, cols = mycolor) + # nolint
    NoLegend() +
    ggtitle("")
  pic2 <- DimPlot(objs, group.by = group_var, split.by = sample_var, label = TRUE, pt.size = 2.0, reduction = reduction, cols = mycolor) + # nolint
    NoLegend() +
    ggtitle("")
  return(list(pic1 = pic1, pic2 = pic2))
}

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
  tbl <- table(objs@meta.data[[group_var]], objs@meta.data[[sample_var]])
  res <- chisq.test(tbl)
  roe <- tbl / res$expected
  pic <- as.data.frame(roe) %>%
    rename(Cell = Var1, sample = Var2, roe = Freq) %>%
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


pdf(paste0(name, "_cellsurvey_outplots.pdf"), width = 26, height = 8)
for (group_var in group_var_list) {
  umap_plots <- create_umap_plots(objs, group_var = group_var, mycolor_string)
  percent_plots <- create_percent_plot(objs, group_var = group_var, mycolor_string)
  roe_plot <- create_roe_plot(objs, group_var = group_var, sample_var = sample_var, mycolor_string = "lightblue,purple")
  if (is.list(umap_plots)) {
    for (plot_name in names(umap_plots)) {
      print(umap_plots[[plot_name]])
    }
  } 
  if (is.list(percent_plots)) {
    for (plot_name in names(percent_plots)) {
      print(percent_plots[[plot_name]])
    }
  }
  print(roe_plot)
}
dev.off()
