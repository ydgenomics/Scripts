# Date: 20250610
# Image: enrich-R--04
# libraries: org.Cthalictroides.eg.db,org.Pcirratum.eg.db,org.Ahypogaea.eg.db
# gene_csv: gene, cluster, p_val_adj

############# input section ###########
# gene_csv="/data/work/0.peanut/orgdb/preprocess.csv"
# kegg_info_RData="/script/build_orgdb/kegg_info.RData"
# db="/data/work/0.peanut/orgdb/output"
# minp=0.05
# genus="Arachis"
# species="hypogaea"
############################################

library(ggplot2)
library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(tidyverse)
library(optparse)

option_list <- list(
  make_option(c("--gene_csv"), type = "character", default = "/data/work/0.peanut/orgdb/preprocess.csv", help = "input the csv of leiden_0.5"),
  make_option(c("--kegg_info_RData"), type = "character", default = "/data/work/0.peanut/orgdb/test5/kegg_info.RData", help = "Kegg info Rdata"),
  make_option(c("--db"),type = "character", default = "/data/work/0.peanut/orgdb/test4",help = "Name of built db for enrich"),
  make_option(c("--minp"), type = "numeric", default = 0.05, help = "filter marker gene limited by min pvalue_adj"),
  make_option(c("--genus"), type = "character", default = "Arachis", help = "Genus name", metavar = "character"),
  make_option(c("--species"), type = "character", default = "hypogaea", help = "Species name", metavar = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
gene_csv <- opt$gene_csv
kegg_info_RData <- opt$kegg_info_RData
db <- opt$db
minp <- opt$minp
genus <- opt$genus
species <- opt$species

# Good for wdl
parent_dir <- db
# library
db_name <- paste0("org.", substr(genus, 1, 1), species, ".eg.db"); print(db_name)
DB <- paste0(parent_dir, '/', db_name); print(DB)
install.packages(DB, repos = NULL, type = "sources")
do.call(library, list(db_name))
db <- get(db_name)
columns(db)


markers <- read.csv(gene_csv, header = TRUE, stringsAsFactors = FALSE)

check_marker_genes <- function(markers, db) {
    required_cols <- c("gene", "cluster", "p_val_adj")
    missing_cols <- setdiff(required_cols, colnames(markers))
    if (length(missing_cols) > 0) {
        stop(paste(
            "Error: The following required columns are missing in gene_csv:",
            paste(missing_cols, collapse = ", ")
        ))
    }
    head(markers) # gene, cluster, p_val_adj

    # Check
    db_gid <- keys(db, keytype = "GID")
    common_genes <- markers$gene[markers$gene %in% db_gid]
    num_common_genes <- length(common_genes)
    total_genes <- length(markers$gene)
    percentage <- (num_common_genes / total_genes) * 100

    cat("First 10 GIDs in database:\n")
    print(head(db_gid, 10))
    cat("Total number of GIDs in database:", length(db_gid), "\n")
    cat("\nFirst 10 genes in input gene_csv:\n")
    print(head(markers$gene, 10))
    cat("Total number of genes in gene_csv:", total_genes, "\n")
    cat("\nNumber of genes present in database:", num_common_genes, "\n")
    cat("Percentage of input genes matched to database:",
            round(percentage, 2), "%\n")
}

check_marker_genes(markers, db)


# pathway and kegg
pathway2gene <- AnnotationDbi::select(db,keys = keys(db),columns = c("Pathway","Ko")) %>%
  na.omit() %>%
  dplyr::select(Pathway, GID)
load(kegg_info_RData)

# Output dictionary
filepath <- paste0(species, "_enrich")
dir.create(filepath)
setwd(filepath)

for(i in unique(markers$cluster)){
    marker_subset <- filter(markers, cluster == i)
    length(marker_subset$gene)
    gene_list <- marker_subset %>% filter(p_val_adj < minp)
    gene_list <- gene_list$gene
    #gene_list <- paste0(gene_list, ".1") # gene_id == GID
    length(gene_list)
    # run enrich
    go_data <- enrichGO(gene = gene_list,OrgDb = db,keyType = 'GID',ont = 'ALL',qvalueCutoff = 0.05,pvalueCutoff = 0.05)
    go_data <- as.data.frame(go_data)
    kegg_result <- enricher(gene_list,TERM2GENE = pathway2gene, TERM2NAME = pathway2name,pvalueCutoff = 0.05,qvalueCutoff = 0.05)
    kegg_data <- as.data.frame(kegg_result); dim(kegg_data)
    if (nrow(go_data) > 0 && nrow(kegg_data) > 0) {
        kegg_data$ONTOLOGY <- "KEGG"
        col_names <- names(kegg_data)
        kegg_data <- kegg_data[, c("ONTOLOGY", col_names[!col_names %in% "ONTOLOGY"])]
        data <- rbind(go_data, kegg_data)
    } else {
        data <- go_data
        print(paste0(i, " lacked enrichment kegg information"))
    }
    if (nrow(data) > 0) {
        print(paste0("Data is not empty, Proceeding ", i))
        length(data$ID)
        data$Name <- paste0(data$ID,"_",data$Description)
        write.table(data, file = paste0(i,"_enrich.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        data_subset <- data %>% arrange(desc(Count)) %>% head(60)
        length(data_subset$ID)
        data_subset <- data_subset %>% mutate(GeneRatio = as.numeric(gsub("/.*", "", GeneRatio)) / as.numeric(gsub(".*/", "", GeneRatio)))
        data_subset$Name <- ifelse(nchar(data_subset$Name) > 100, substr(data_subset$Name, 1, 100), data_subset$Name)
        # viusal1
        if (length(data_subset$ID) > 0) {
        pdf(paste0(i,"_plot1.pdf"))
        plot1 <- ggplot(data_subset, aes(y = GeneRatio, x = reorder(Name, GeneRatio))) + 
            geom_bar(stat = "identity", aes(fill = p.adjust), width = 0.8) +  
            scale_fill_gradient(low = "red", high = "blue") +  
            facet_grid(ONTOLOGY ~ ., scales = "free", space = "free") +  
            coord_flip() + xlab("Name") + ylab("GeneRatio") + labs(title = paste0("Group ", i, " GO and KEGG Enrich")) + 
            theme(
                axis.text.x = element_text(size = 10), 
                axis.text.y = element_text(size = 5), 
                axis.title.x = element_text(size = 12),  
                axis.title.y = element_text(size = 12)) +
            geom_text(aes(label = Count), vjust = 0, size = 1.5) +
            scale_size_continuous(range = c(0.1, 3)) 
        print(plot1)
        dev.off()}
    } else {
        print("Data is empty. Skipping the code.")
    }
}