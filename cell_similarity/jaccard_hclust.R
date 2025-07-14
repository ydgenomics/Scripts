# Date: 20250714
# Image: metaNeighbor-R--02 /opt/conda/bin/R
# 做了jaccard相似性的计算基于SCT的FindAllMarkers的基因和hclust基于SCT的counts的基因表达
# Â© EMBL-European Bioinformatics Institute, 2023
# Yuyao Song <ysong@ebi.ac.uk>

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggdendro)
library(ape)
library(glmGamPoi)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input_file"),
    type = "character", default = "/data/work/script/rename0626/D1/D1_0626.rds",
    help = "Path to input file"
  ),
  make_option(c("-o", "--output_name"),
    type = "character", default = "D1",
    help = "Output file prefix name"
  ),
  make_option(c("-b", "--batch_key"),
    type = "character", default = "sample",
    help = "Batch key for integration"
  ),
  make_option(c("-c", "--cluster_key"),
    type = "character", default = "rename0626",
    help = "Cluster key for integration"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))
input_file <- opt$input_file
out_put_name <- opt$output_name
batch_key <- opt$batch_key
cluster_key <- opt$cluster_key

# input_file <- "/data/work/script/rename0626/D1/D1_0626.rds" 
# output_name <- "D1" 
# batch_key <- "sample" 
# cluster_key <- "rename0626"

scRNA <- readRDS(input_file)
scRNA
colnames(scRNA@meta.data)
#scRNA <- subset(scRNA, cells = sample(Cells(scRNA), 500))
# Preprocess: Check whether existing SCT
if (!"SCT" %in% Assays(scRNA)) {
    print("SCT layer not existing, so run SCTransform and integrateLayers.")
    scRNA[["RNA"]] <- split(scRNA[["RNA"]], f = scRNA@meta.data[[batch_key]]) # this command will split object in layers!
    scRNA <- SCTransform(scRNA, vst.flavor = "v2")
    scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
    scRNA <- IntegrateLayers(object = scRNA, method = HarmonyIntegration,
                             orig.reduction = "pca", new.reduction = 'harmony',
                             assay = "SCT", verbose = FALSE)
}

DefaultAssay(scRNA) <-"SCT"
Idents(scRNA) <- paste0(scRNA@meta.data[[cluster_key]],"_of_",scRNA@meta.data[[batch_key]])

# Find Idents() All markers
scRNA <- PrepSCTFindMarkers(scRNA)
df_gene=FindAllMarkers(scRNA,only.pos = T,logfc.threshold = 0.1)
table(df_gene$cluster)

# Calculate
cluster=names(table(df_gene$cluster))
## jaccard
df_ja=c()
for (i in cluster) {
  ja=c()
  for (j in cluster) {
    a=df_gene[df_gene$cluster==i,]
    a=a$gene
    b=df_gene[df_gene$cluster==j,]
    b=b$gene
    jaccard=length(intersect(a,b))/length(union(a,b))
    ja=c(ja,jaccard)
  }
  df_ja=rbind(df_ja,ja)
}
rownames(df_ja)=cluster
colnames(df_ja)=cluster
head(df_ja)
# df_ja2
df_ja2 <- as.data.frame(df_ja)
df_ja2$cluster1 <- rownames(df_ja2)
df_ja2 <- tidyr::pivot_longer(df_ja2,!cluster1, names_to="cluster2",values_to ="jaccard")
head(df_ja2)

# this is cribbed from 
# https://stackoverflow.com/questions/42047896/joining-a-dendrogram-and-a-heatmap
# to align dendrogram with dotplot
#scale_rows <- function (x) {
#  m = apply(x, 1, mean, na.rm = T)
#  s = apply(x, 1, sd, na.rm = T)
#  return((x - m)/s)
#}
#mat = switch('row', none = df_ja, row = scale_rows(df_ja), column = t(scale_rows(t(df_ja))))
d = dist(df_ja, method = 'euclidean') #计算距离矩阵
ddgram = hclust(d, method = 'complete') #层次聚类

#ddgram <- hclust(dist(df_ja))
ddata <- dendro_data(ddgram, type = 'rectangle') # extract into lists of data
gene_pos_table <- with(ddata$labels, data.frame(y_center = x, gene = as.character(label), height = 1))
# axis munging <- This is where the magic happens
gene_axis_limits <- with(
  gene_pos_table, 
  c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))) +  0.1 * c(-1, 1)

ddata <- with(
  segment(ddata), 
  data.frame(x = y, y = x, xend = yend, yend = xend))

fancy_tree_plot <-  ggplot((ddata)) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  scale_x_reverse(expand = c(0, 0.5)) + 
  scale_y_continuous(breaks = gene_pos_table$y_center, 
                     labels = gene_pos_table$gene, 
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "Distance", y = "", colour = "", size = "") +
  theme_dendro() 


dotplot <- df_ja2 %>%
  mutate(cluster1 = factor(cluster1, levels = gene_pos_table$gene)) %>% 
  mutate(cluster2 = factor(cluster2, levels = gene_pos_table$gene)) %>% 
  ggplot(aes(x=cluster1, y = cluster2, color = jaccard, size = jaccard)) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
#  theme(axis.line  = element_blank()) +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme_bw() +
  scale_color_distiller(
    palette = 'Reds',
    direction = 1,
    name = 'Log-normalised\nexpression',
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"),limits = c(0,0.5),oob = scales::squish
  )+
#  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,0.5),oob = scales::squish)+
  scale_y_discrete(position = "right")+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#################################################

p <- plot_grid(fancy_tree_plot, NULL, dotplot, nrow = 1, rel_widths = c(0.5,-0.1, 2), align = 'h')
p
#ggsave2(p,file=paste0("Jaccard_",out_put_name,".pdf"),width=20,height=15)
file_name <- paste("Jaccard_", out_put_name, ".pdf", sep = "")
number <- length(unique(Idents(scRNA)))/2
ggsave2(p, file = file_name, width = number, height = number)
sample.cluster <- AggregateExpression(scRNA,group.by = "ident",slot="counts")$SCT
#sample.cluster <- AggregateExpression(scRNA, group.by = "ident", layer="counts")$SCT
hc = hclust(dist(t(as.matrix(sample.cluster))))
#save(sample.cluster,hc,file="hc.RData")
#load(file="hc.RData")
dend = as.dendrogram(hc)

#.libPaths(c("~/R/x86_64-conda-linux-gnu-library/4.1"))
library(dendextend)
library(circlize)

clusM <- c(sapply(strsplit(as.character(unique(Idents(scRNA))),'_of_'), "[", 2))
names(clusM)<-unique(Idents(scRNA))

colors <- c("lightcoral","lightseagreen","grey","red")
names(colors) <- unique(clusM)

cols <- c()
sname <- c()
for (name in labels(dend))
   {sname <- c(sname,name)
    tissue <- clusM[which(names(clusM)==name)] 
    color <- colors[which(names(colors)==tissue)]
    cols <- c(cols,color)}
file_name2 <- paste("hclust_", out_put_name, ".pdf", sep = "")
pdf(file_name2,height=number,width=number)
#plot(as.phylo(hc), type = "fan",tip.color = cols,
#     label.offset = 1, cex = 0.7)

dend <- dend %>% set("labels_col", cols) %>% # change color
  set("labels_cex", 0.5) %>% # Change size
  set("branches_lwd", 2) %>%
  set("branches_k_color", k = 4)
plot(dend) # plot
ggd1 <- as.ggdend(dend)
#number_of_bar <- length(labels(dend))
#angle <-  90 - 360 * (c(1:length(labels(dend)))-0.5) /number_of_bar
#hjust<-ifelse( angle < -90, 1, 0)
#angle<-ifelse(angle < -90, angle+180, angle)
circlize_dendrogram(dend,labels_track_height =0.5,dend_track_height = 0.3) 
#ggplot(ggd1) + 
# scale_y_reverse(expand = c(0.2, 0)) +
#  coord_polar(theta="x")+
#  theme(axis.text.x = element_text(
#    angle= -90 - 360 / length(labels(dend)) * seq_along(labels(dend))))
#  coord_radial(rotate_angle = TRUE, expand = FALSE)
dev.off()

#obj <- JoinLayers(obj)
#obj [["RNA"]] <- JoinLayers(obj [["RNA"]])

# Assay RNA changing from Assay5 to Assay
#obj[["RNA"]] <- as(obj[["RNA"]], "Assay")

saveRDS(scRNA, file = paste0(out_put_name, "_sct.rds"))