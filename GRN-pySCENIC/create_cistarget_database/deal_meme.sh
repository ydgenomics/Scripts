# Title: deal_genome.R
# Date: 20250522
# Coder: ydgenomics
# Description:
# Input: gtf file and fasta file
# Output:
# Image: 
# Reference: https://mp.weixin.qq.com/s/7-vKrLiFS4Tlkt-rHxEGeQ; https://github.com/aertslab/create_cisTarget_databases

cd /data/work/0.peanut/GRN/output
# 1.使用R代码在PlantTFDB下载合适meme文件
less -S /data/work/0.peanut/GRN/input/Ath_TF_binding_motifs.meme
# 2.对.meme进行处理，得到每个 motif 对应一个矩阵文件，需要以.cb 结尾，文件名和 motif 保持一致
## 提取motif及矩阵
grep -E "MOTIF|^[[:space:]]*[0-9]" /data/work/0.peanut/GRN/input/Ath_TF_binding_motifs.meme | sed 's/MOTIF />/g' | sed 's/^[[:space:]]*//g' > tf_motif_matrix.txt
## 替换空格
sed -i 's/>MA\([0-9.]*\) />MA\1_/' tf_motif_matrix.txt
less tf_motif_matrix.txt
## 保证名字为motif1的名字
awk '/^>/ {print ">" $2; next} {print}' tf_motif_matrix.txt > formatted_tf_motif_matrix.txt
less formatted_tf_motif_matrix.txt
# 3.输出保存到每个文件
mkdir motif_dir
cd motif_dir
awk '/^>/{if(file) close(file); filename=substr($0,2)".cb"; print $0 > filename; file=filename; next} {print >> file}'  ../formatted_tf_motif_matrix.txt
## 准备每个motif id文件
grep ">" /data/work/0.peanut/GRN/output/formatted_tf_motif_matrix.txt|sed 's/>//g' > /data/work/0.peanut/GRN/output/motifs_id.txt
head /data/work/0.peanut/GRN/output/motifs_id.txt