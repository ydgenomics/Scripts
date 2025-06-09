
version 1.0
workflow Enrich_scRNAseq_V1_1_0_4{
  input{
    File emapper_xlsx
    File gene_csv
    File ic_tsv
    File relations_full_tsv
    File go_obo
    File gofigure_py
    String genus
    String species
    String taxid
    Float minp=0.05
  }
  String url_enrich = "public-library/yangdong_cbd986609f064080b930131023416378_public:latest" # SoupX-R--01
  String url_gofigure = "stereonote_hpc/yangdong_62c01e5f5e724b32b2591c91a9a722e6_private:latest"
  call build_orgdb{
    input:
    emapper_xlsx=emapper_xlsx,
    genus=genus,
    species=species,
    taxid=taxid,
    cpu=2,
    mem=8,
    url=url_enrich,
  }
  call check{
    input:
    DB=build_orgdb.db,
    gene_csv=gene_csv,
    genus=genus,
    species=species,
    cpu=1,
    mem=4,
    url=url_enrich,
  }
  call enrich{
    input:
    DB=build_orgdb.db,
    gene_csv=gene_csv,
    genus=genus,
    species=species,
    minp=minp,
    cpu=2,
    mem=8,
    url=url_enrich,
  }
  call gofigure{
    input:
    result_txt=enrich.result_txt,
    ic_tsv=ic_tsv,
    relations_full_tsv=relations_full_tsv,
    go_obo=go_obo,
    gofigure_py=gofigure_py,
    max_label=15,
    cpu=2,
    mem=8,
    url=url_gofigure,
  }
  output{
    File result=build_orgdb.db
    File csv=check.csv
    File pdf=enrich.result
    Array[File] enrich_txt=enrich.result_txt
    File gofigure_result=gofigure.result
  }
}
task build_orgdb{
  input {
    File emapper_xlsx
    String genus
    String species
    String taxid
    Int cpu
    Int mem
    String url
  }
  command <<<
    mkdir result
    cd result
    /opt/conda/bin/Rscript /script/build_orgdb/build_orgdb.R \
    --emapper_xlsx ~{emapper_xlsx} --taxid ~{taxid} --genus ~{genus} --species ~{species}
  >>>
  runtime {
    docker_url: "~{url}"
    req_cpu: cpu
    req_memory: "~{mem}Gi"
  }
  output {
    File db = "result"
  }
}

task check{
  input {
    File DB
    File gene_csv
    String genus
    String species
    Int cpu
    Int mem
    String url
  }
  command <<<
    /opt/conda/bin/R --vanilla --slave <<EOF
    genus <- "~{genus}"
    species <- "~{species}"
    gene_csv <- "~{gene_csv}"
    DB <- "~{DB}"
    # gene_set
    gene_set <- read.csv(gene_csv, header = TRUE, sep = ",")
    head(gene_set)
    print(head(gene_set$gene, 10))
    ############################ make gene map to gid ##############################
    #gene_set$gene <- paste0(gene_set$gene, ".1") # gene_id == GID
    ################################################################################
    write.csv(gene_set, "preprocess.csv")
    EOF
  >>>
  runtime {
    docker_url: "~{url}"
    req_cpu: cpu
    req_memory: "~{mem}Gi"
  }
  output {
    File csv="preprocess.csv"
  }
}

task enrich{
  input {
    File DB
    File gene_csv
    String genus
    String species
    Float minp
    Int cpu
    Int mem
    String url
  }
  command <<<
    mkdir enrich
    cd enrich
    /opt/conda/bin/R --vanilla --slave <<EOF
    # Date: 20250609
    # Image: enrich-R--04
    # libraries: org.Cthalictroides.eg.db,org.Pcirratum.eg.db,org.Ahypogaea.eg.db
    # gene_csv: gene, cluster, p_val_adj

    library(ggplot2)
    library(tidyverse)
    library(dplyr)
    library(clusterProfiler)
    library(optparse)

    option_list <- list(
    make_option(c("--gene_csv"), type = "character", default = "~{gene_csv}", help = "input the csv of leiden_0.5"),
    make_option(c("--minp"), type = "numeric", default = ~{minp}, help = "filter marker gene limited by min pvalue_adj"),
    make_option(c("--db"),type = "character", default = "~{DB}",help = "Name of built db for enrich"),
    make_option(c("--genus"), type = "character", default = "~{genus}", help = "Genus name", metavar = "character"),
    make_option(c("--species"), type = "character", default = "~{species}", help = "Species name", metavar = "character")
    )
    opt <- parse_args(OptionParser(option_list = option_list))

    # Good for wdl
    parent_dir <- opt$db
    paths <- list.files(parent_dir, full.names = TRUE, recursive = FALSE)
    print(paths)
    DB <- paths[1]

    # library
    db_name <- paste0("org.", substr(opt$genus, 1, 1), opt$species, ".eg.db")
    print(db_name)
    install.packages(DB, repos = NULL, type = "sources")
    do.call(library, list(db_name))
    db <- get(db_name)
    columns(db)


    markers <- read.csv(opt$gene_csv, header = TRUE, stringsAsFactors = FALSE)
    head(markers) # gene, cluster, p_val_adj

    check_marker_genes <- function(markers, db) {
        required_cols <- c("gene", "cluster", "p_val_adj")
        missing_cols <- setdiff(required_cols, colnames(markers))
        if (length(missing_cols) > 0) {
            stop(paste(
                "Error: The following required columns are missing in gene_csv:",
                paste(missing_cols, collapse = ", ")
            ))
        }
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
    load("/data/work/0.peanut/orgdb/kegg_info.RData")

    # Output dictionary
    filepath <- paste0(opt$species, "_enrich")
    dir.create(filepath)
    setwd(filepath)

    for(i in unique(markers$cluster)){
        marker_subset <- filter(markers, cluster == i)
        length(marker_subset$gene)
        gene_list <- marker_subset %>% filter(p_val_adj < opt$minp)
        gene_list <- gene_list$gene
        #gene_list <- paste0(gene_list, ".1") # gene_id == GID
        length(gene_list)
        # run enrich
        go_data <- enrichGO(gene = gene_list,OrgDb = db,keyType = 'GID',ont = 'ALL',qvalueCutoff = 0.05,pvalueCutoff = 0.05)
        go_data <- as.data.frame(go_data)
        kegg_result <- enricher(gene_list,TERM2GENE = pathway2gene, TERM2NAME = pathway2name,pvalueCutoff = 0.05,qvalueCutoff = 0.05)
        kegg_data <- as.data.frame(kegg_result)
        dim(kegg_data)
        kegg_data$ONTOLOGY <- "KEGG"
        col_names <- names(kegg_data)
        kegg_data <- kegg_data[, c("ONTOLOGY", col_names[!col_names %in% "ONTOLOGY"])]
        if (nrow(go_data) > 0 && nrow(kegg_data) > 0) {
            data <- rbind(go_data, kegg_data)
        } else {
            data <- go_data
            print(paste0(i, " lacked enrichment kegg information"))
        }
        if (nrow(data) > 0) {
            print(paste0("Proceeding: ", i))
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
            print(paste0("Data of", i, "is empty. Skipping the code."))
        }
    }
    EOF
  >>>
  runtime {
    docker_url: "~{url}"
    req_cpu: cpu
    req_memory: "~{mem}Gi"
  }
  output {
    File result="enrich"
    Array[File] result_txt = glob("enrich/*.txt")
  }
}

task gofigure{
  input {
    Array[File] result_txt
    File ic_tsv
    File relations_full_tsv
    File go_obo
    File gofigure_py
    Int max_label
    Int cpu
    Int mem
    String url
  }
  command <<<
    for c in ~{sep="," result_txt}; do
        echo $c >> result_txt.txt
        basename=$(basename "$c" _enrich.txt)
        echo $basename >> result_name.txt
    done
    /opt/software/R/bin/R --vanilla --slave <<EOF
    file_path <- "result_txt.txt"
    file_content <- readLines(file_path, warn = FALSE)
    file_list <- strsplit(file_content, ",")[[1]]
    file_list <- trimws(file_list)
    print(file_list)

    output_files <- c()
    for (file_path in file_list) {
    # 提取文件名
    file_name <- basename(file_path)
    
    # 提取文件名的前缀（不包括后缀）
    prefix <- tools::file_path_sans_ext(file_name)
    
    # 构建输出文件名
    output_file_path <- paste0(prefix, "_output_standard_gofigure_input.tsv")
    
    # 读取文件内容
    result <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)
    
    # 提取 ID 和 p.adjust 列
    result <- result[, c("ID", "p.adjust")]
    
    # 筛选 ID 列中以 "GO:" 开头的行
    result <- result[grep("^GO:", result$ID), ]
    
    # 重命名列
    colnames(result) <- c("GOterm", "enrichment_P-value")
    
    # 将每一行的内容合并为一个字符串，列之间用制表符分隔
    result_lines <- apply(result, 1, function(row) {
        paste(row, collapse = "\t")
    })
    
    # 将所有行的内容合并为一个字符串，行之间用换行符分隔
    result_text <- paste(result_lines, collapse = "\n")
    
    # 保存到文件
    write(result_text, file = output_file_path)
    
    # 将输出文件路径转换为绝对路径
    absolute_output_file_path <- normalizePath(output_file_path)
    
    # 打印输出文件路径
    print(paste("Output saved to:", absolute_output_file_path))
    
    # 将绝对路径保存到 output_files 向量中
    output_files <- c(output_files, absolute_output_file_path)
    }
    # 将所有输出文件名保存到一个文件中
    output_summary_file <- file.path("output_standard_gofigure_input.txt")
    write(paste(output_files, collapse = ","), file = output_summary_file)
    EOF

    mkdir gofigure_result
    cd gofigure_result

    result_name_file="../result_name.txt"
    output_file="../output_standard_gofigure_input.txt"

    # 读取文件内容并按逗号分割
    IFS=',' read -r -a names <<< "$(cat "$result_name_file")"
    IFS=',' read -r -a outputs <<< "$(cat "$output_file")"

    # 获取数组长度
    len_names=${#names[@]}
    len_outputs=${#outputs[@]}

    # 确保两个数组长度一致
    if [ "$len_names" -ne "$len_outputs" ]; then
    echo "Error: The number of names and outputs does not match."
    exit 1
    fi

    # 循环处理每个名称和对应的输出文件
    i=0
    while [ $i -lt $len_names ]; do
    name=${names[$i]}
    output=${outputs[$i]}

    echo "Processing $name with output file $output"
    mkdir "$name"
    tsv_path="$output"
    max_label=~{max_label}

    # run go-figure
    mkdir "./$name/data"
    cp ~{ic_tsv} "$name/data"
    cp ~{relations_full_tsv} "$name/data"
    cp ~{go_obo} "$name/data"
    mkdir "$name/go-figure-output"
    /software/miniconda/envs/go-figure/bin/python ~{gofigure_py} \
        -i "$tsv_path" -j standard -m "$max_label" -o "$name/go-figure-output"
    rm -r "$name/data"

    # 自增索引
    i=$((i + 1))
    done
  >>>
  runtime {
    docker_url: "~{url}"
    req_cpu: cpu
    req_memory: "~{mem}Gi"
  }
  output {
    File result="gofigure_result"
  }
}