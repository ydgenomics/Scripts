# Title: jaccard_hclust.R
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
    type = "character", default = NULL,
    help = "Path to input file"
  ),
  make_option(c("-o", "--output_name"),
    type = "character", default = NULL,
    help = "Output file prefix name"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input_file)){
  opt$input_file <- "/data/work/integration/input/Peanut-unsoupx.cg.rds"
} 
if (is.null(opt$output_name)){
  opt$output_name <- "data"
} 
input_file <- opt$input_file
out_put_name <- opt$output_name
batch_key <- "biosample"
cluster_key <- "leiden_res_0.50"


scRNA=readRDS(input_file)
scRNA
#scRNA <- subset(scRNA, cells = sample(Cells(scRNA), 500))
# 如果输入对象不包含SCT层才做这个处理
if (!"SCT" %in% Assays(scRNA)) {
    print("SCT layer not existing, so run SCTransform and integrateLayers.")
    scRNA[["RNA"]] <- split(scRNA[["RNA"]], f = scRNA@meta.data[batch_key])
    scRNA <- SCTransform(scRNA, vst.flavor = "v2")
    scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
    scRNA <- IntegrateLayers(object = scRNA, method = HarmonyIntegration,
                      orig.reduction = "pca", new.reduction = 'harmony',
                       assay = "SCT", verbose = FALSE)
}


DefaultAssay(object = scRNA) <-"SCT"
#don't found scRNA$celltype, so use leiden_res_0.50 replaced celltype
#### in order to solve the empty problem of celltype
if ("celltype" %in% colnames(scRNA@meta.data)) {
  if (all(is.na(scRNA@meta.data$celltype))) {
    scRNA@meta.data$celltype <- scRNA@meta.data$leiden_res_0.50
  }
} else {
  scRNA@meta.data$celltype <- scRNA@meta.data$leiden_res_0.50
}

Idents(scRNA) <- paste0(scRNA$celltype,"_of_",scRNA$batch)

scRNA <- PrepSCTFindMarkers(scRNA)
df_gene=FindAllMarkers(scRNA,only.pos = T,logfc.threshold = 0.1)
#table(df_gene$cluster)

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
df_ja2 <- as.data.frame(df_ja)
df_ja2$cluster1 <- rownames(df_ja2)
df_ja2 <- tidyr::pivot_longer(df_ja2,!cluster1, names_to="cluster2",values_to ="jaccard")

# this is cribbed from 
# https://stackoverflow.com/questions/42047896/joining-a-dendrogram-and-a-heatmap
# to align dendrogram with dotplot
#scale_rows <- function (x) {
#  m = apply(x, 1, mean, na.rm = T)
#  s = apply(x, 1, sd, na.rm = T)
#  return((x - m)/s)
#}
#mat = switch('row', none = df_ja, row = scale_rows(df_ja), column = t(scale_rows(t(df_ja))))
d = dist(df_ja, method = 'euclidean')
ddgram = hclust(d, method = 'complete')

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
ggsave2(p, file = file_name, width = 20, height = 15)
sample.cluster<-AggregateExpression(scRNA,group.by = "ident",slot="counts")$SCT
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
pdf(file_name2,height=8,width=8)
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