library(rtracklayer)
library(tidyverse)
library(Biostrings)

gtf <- import("/data/input/Files/husasa/Ref/arahy.Tifrunner.gnm2.ann2.PVFB.gene_models_main.gtf", format = "gtf") %>%
  as.data.frame()

head(gtf)
#> colnames(gtf)
# [1] "seqnames"      "start"         "end"           "width"        
# [5] "strand"        "source"        "type"          "score"        
# [9] "phase"         "transcript_id" "gene_id"
#> unique(gtf$type)
#[1] transcript exon       CDS

# get gene promoters regions
genes <- gtf %>%
  dplyr::filter(type == "CDS" & gene_id != "NA") %>%
  dplyr::select(seqnames, start, end, strand, gene_id) %>%
  dplyr::mutate(start2 = ifelse(strand == "+", start - 3000, end + 1),
                end2 = ifelse(strand == "+", start - 1, end + 3000))
#对于正链（strand == "+"），启动子区域从 start - 3000 到 start - 1。
#对于负链（strand == "-"），启动子区域从 end + 1 到 end + 3000。
genes$seqnames <- as.character(genes$seqnames)
genes$seqnames <- sapply(strsplit(genes$seqnames, split = "\\."), "[", 4)
head(genes$seqnames, 5)
head(genes) # genes$seqnames必须要和names(genome)一致

########## Rename chromosomes of genome ##########
# load genome
genome <- readDNAStringSet("/data/input/Files/husasa/Ref/arahy.Tifrunner.gnm2.J5K5.genome_main.fa")
genome

# rename seqnames
#names(genome) <- sapply(strsplit(names(genome),split = " "), "[", 1) #unmatch peanut genome
names(genome) <- sapply(strsplit(names(genome), split = "\\."), "[", 4)
head(names(genome), 25)

########### Cycle extracting sequences of promoter ###########
# extract promter sequences
# x = 1
lapply(1:nrow(genes), function(x) {
  print(x)
  tmp <- genes[x, ]
  # try extract seq
  out <- tryCatch(
    seq <- genome[[tmp$seqnames]][tmp$start2:tmp$end2],
    error = function(e){ return(NULL) }
  )
  # check
  if(!is.null(out)){
    # strand
    if(tmp$strand == "-"){
      out <- reverseComplement(out)
    }else{
      out
    }
  }else{
    DNAString()
  }
}) %>% DNAStringSet() -> seq_list

names(seq_list) <- make.unique(genes$gene_id)

seq_list

# filte seq length
fasta_filtered <- seq_list[width(seq_list) > 1]
head(fasta_filtered)
# output fasta format
writeXStringSet(fasta_filtered, filepath = "3kpromoter.fasta", format = "fasta")