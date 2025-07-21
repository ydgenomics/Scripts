[同源基因比对：顶刊都在用的跨物种分析](https://mp.weixin.qq.com/s/jv2Z8NVWZwzeVjmm5c9NNg)
[使用ensemble得到两个物种的同源基因：从BioMart中下载拟南芥和小麦的直系同源基因（Orthologues）](https://www.jianshu.com/p/5de2c98797f2)

同源基因和直系同源基因的区别？

学习使用blast做两个不同物种的蛋白质序列找同源基因
BLASTP本身主要用于单序列与数据库的比对，而不是用于多序列比对
[史上最详细的blast安装附视频](https://mp.weixin.qq.com/s/rEBqjN-fGOp_loTmyEuMJA)
```shell
########## Install blastp ##########
cd /software
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz
tar -zxvf ncbi-blast-2.16.0+-x64-linux.tar.gz
# environment location: /software/ncbi-blast-2.16.0+/bin
vim ~/.bashrc
export PATH=/software/ncbi-blast-2.16.0+/bin:$PATH
source ~/.bashrc
blastp -h
rm ncbi-blast-2.16.0+-x64-linux.tar.gz
```
```shell
########## Using blastp ##########
subject_fasta="/data/work/0.peanut/GRN/input/Ath_pep.fas"
query_fasta="/data/work/0.peanut/GRN/input/arahy.Tifrunner.gnm2.ann2.PVFB.protein.faa"
#cd /data/work/0.peanut/GRN/input/arabidopsis_db
# Make database of AT/Others
makeblastdb -in $subject_fasta -dbtype prot -out arabidopsis_db
# Query input fasta
blastp -query $query_fasta -db arabidopsis_db -out blastp_results.txt -outfmt 6 -evalue 1e-5
```

blastp不适合多序列比对，使用Clustal加快序列比对
[Clustal Omega—广泛使用的多序列比对工具](https://mp.weixin.qq.com/s/f9pEFWJJoNCqlFEfd77aOA)

```shell
blastp -help
USAGE
  blastp [-h] [-help] [-import_search_strategy filename]
    [-export_search_strategy filename] [-task task_name] [-db database_name]
    [-dbsize num_letters] [-gilist filename] [-seqidlist filename]
    [-negative_gilist filename] [-negative_seqidlist filename]
    [-taxids taxids] [-negative_taxids taxids] [-taxidlist filename]
    [-negative_taxidlist filename] [-no_taxid_expansion] [-ipglist filename]
    [-negative_ipglist filename] [-entrez_query entrez_query]
    [-db_soft_mask filtering_algorithm] [-db_hard_mask filtering_algorithm]
    [-subject subject_input_file] [-subject_loc range] [-query input_file]
    [-out output_file] [-evalue evalue] [-word_size int_value]
    [-gapopen open_penalty] [-gapextend extend_penalty]
    [-qcov_hsp_perc float_value] [-max_hsps int_value]
    [-xdrop_ungap float_value] [-xdrop_gap float_value]
    [-xdrop_gap_final float_value] [-searchsp int_value] [-seg SEG_options]
    [-soft_masking soft_masking] [-matrix matrix_name]
    [-threshold float_value] [-culling_limit int_value]
    [-best_hit_overhang float_value] [-best_hit_score_edge float_value]
    [-subject_besthit] [-window_size int_value] [-lcase_masking]
    [-query_loc range] [-parse_deflines] [-outfmt format] [-show_gis]
    [-num_descriptions int_value] [-num_alignments int_value]
    [-line_length line_length] [-html] [-sorthits sort_hits]
    [-sorthsps sort_hsps] [-max_target_seqs num_sequences]
    [-num_threads int_value] [-mt_mode int_value] [-ungapped] [-remote]
    [-comp_based_stats compo] [-use_sw_tback] [-version]
num_threads 
```