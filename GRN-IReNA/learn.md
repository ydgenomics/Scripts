[树棉：*Gossypium arboreum*](https://baike.baidu.com/item/%E6%A0%91%E6%A3%89/1706952?fromModule=search-result_lemma)

计算树棉基因组大小`awk '{sum += length($0)} END {print sum}' /data/input/Files/taoziyi/cotton_atac/NB2025053011270768166314/genome/genome.fa` **1444625381**

[*Gossypium arboreum* of TF_blinding_motif](https://planttfdb.gao-lab.org/download/motif/Gar_TF_binding_motifs_information.txt)
```R
download.file(
  url = "https://planttfdb.gao-lab.org/download/motif/Gar_TF_binding_motifs_information.txt",
  destfile = "Gar_TF_binding_motifs_information.txt",
  mode = "wb"
)
```