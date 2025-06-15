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

构建cotton的`motif1`
```R
> head(motif1)
  Accession      ID   Name         TFs                             EnsemblID
1    M00001 MYOD_01   MyoD       Myod1                    ENSMUSG00000009471
2    M00002  E47_01    E47        Tcf3                    ENSMUSG00000020167
3    M00004 CMYB_01  c-Myb         Myb                    ENSMUSG00000019982
4    M00005  AP4_01   AP-4       Tfap4                    ENSMUSG00000005718
5    M00006 MEF2_01 MEF-2A Mef2a;Mef2c ENSMUSG00000030557;ENSMUSG00000005583
6    M00007 ELK1_01  Elk-1        Elk1                    ENSMUSG00000009406
```
这个是motif_blinding_tf的信息，最后一列即5列为`gene`
```R
> head(rownames(seurat_object), n=10)
 [1] "Ga14g01907" "Ga05g03180" "Ga13g00005" "Ga14g01554" "Ga14g01555"
 [6] "Ga14g01813" "Ga14g00262" "Ga01g00010" "Ga01g00022" "Ga01g00004"
 ```

 E1_to_Ga的query基因为什么还有不同的后缀，有什么特殊意义吗？？