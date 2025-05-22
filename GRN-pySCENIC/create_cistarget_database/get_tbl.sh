# 1.从PlantTFDB下载对应的TF_binding_motifs_information.txt文件
# 2.数据结构处理输出.tbl文件
awk '
BEGIN {
    # 打印表头
    print "#motif_id\tmotif_name\tmotif_description\tsource_name\tsource_version\tgene_name\tmotif_similarity_qvalue\tsimilar_motif_id\tsimilar_motif_description\torthologous_identity\torthologous_gene_name\torthologous_species\tdescription"
}
NR > 1 {  # 跳过表头行
    # 提取字段
    gene_id = $1
    matrix_id = $3
    # 打印目标格式
    print matrix_id "\t" matrix_id "\t" gene_id "\tPlantTFDB\t5.0\t" gene_id "\t0\tNone\tNone\t1\tNone\tNone\tgene is directly annotated"
}' /data/work/0.peanut/GRN/input/Ath_TF_binding_motifs_information.txt > /data/work/0.peanut/GRN/output/Ath_TF_binding_motifs_information.tbl