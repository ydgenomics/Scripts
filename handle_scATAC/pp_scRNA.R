# Ref: https://jiang-junyao.github.io/IReNA/scRNA-seq-preprocessing

### Read seurat_object
library(IReNA)
seurat_object <- readRDS('seurat_object.rds')
### Calculate the pseudotime and return monocle object
monocle_object <- get_pseudotime(seurat_object,gene.use = rownames(seurat_object))
###Add pseudotime to the Seurat object
### This function only support monocle object from monocle2
seurat_with_time <- add_pseudotime(seurat_object, monocle_object)

### Identify DEGs across pseudotime (qvalue < 0.05 and num_cells_expressed > 0.1)
library(monocle)
monocle_object <- detectGenes(monocle_object, min_expr = 1)
monocle_object <- estimateDispersions(monocle_object)
diff1 <- monocle::differentialGeneTest(monocle_object,fullModelFormulaStr = "~Pseudotime",relative_expr = TRUE)
sig_genes <- subset(diff1, qval < 0.05)
sig_genes <- subset(sig_genes, num_cells_expressed > 0.1)
### select Candidate TFs.
Candidate_TFs <- c()
for (i in 1:nrow(motif1)) {
  gene1 <- strsplit(motif1[i,5],';')[[1]]
  Candidate_TFs <- c(Candidate_TFs,gene1)
}
### Identify expressed TFs
### Canidate TFs in our motif database are ensemble ID, If your gene names in 
### seurat object are Symbol ID, you need to transfer 
### Canidate TFs to Symbol ID first.
expressed_tf <- rownames(extract_expressed_TFs(seurat_object,Candidate_TFs))
expressed_tf <- expressed_tf[!expressed_tf%in%rownames(sig_genes)]
### Refine the seurat object
refined_seurat <- subset(seurat_with_time, features = c(expressed_tf,rownames(sig_genes)))

# ### DEGs used here is the character class
# seurat_with_time <- subset(seurat_with_time, features = DEGs) motif1 