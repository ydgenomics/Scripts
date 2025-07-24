# Enrichment analysis and Visualization(Enrich)
- **Brief:** Bioinformatics provides statistical information such as transcript numbers, but understanding the biological significance behind these numbers is crucial. Enrichment analysis helps reveal biological meanings, such as GO terms and KEGG pathways. Since biological information is highly specialized, effective visualization is also essential.
- **Log:**
  - 1.0.4
    - 250724 支持循环运行多个marker基因列表
  - 20250711 Added the check query column of eggmapper, deleting repeat rows.
  - 20250616 update wdl.md(description)
  - Visualization with enrichplot will be optimized later; GSEA analysis will be added.
  - For building orgDb, input should be `go_obo` for consistency and authority, as required by Go-Figure visualization.
  - The `check` step ensures input gene sets match the background database. Since they may not match exactly, customization (e.g., adding suffixes) may be needed. The next `task: enrich` requires a CSV with columns: `gene_id`, `cluster`, `p_val_adj`. Adjust the `check` step to ensure this structure.
- **Tradition:** enrich_scRNAseq
---

# Input
- **Parameters:**  
*If your work place isn't 035 project, you need download these files form website and submit them to cloud platform*
  - `ko_json` [download](https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001&format=json&filedir=)
  - `go_obo` [download](https://gitlab.com/evogenlab/GO-Figure/-/tree/master/data?ref_type=heads)
  - `ic_tsv` [download](https://gitlab.com/evogenlab/GO-Figure/-/blob/master/data/ic.tsv?ref_type=heads)
  - `relations_full_tsv` [download](https://gitlab.com/evogenlab/GO-Figure/-/commit/48180848bace51314e5dcc4819cc6eb08ab92e45)
---

# Output

---

# Workflow
- **Overview:**  
  1. build orgdb
  2. enrich including Go and KEGG
  3. Visual with go-figure
  4. Summary
   > Note: 不能直接将建好的orgdb作为流程的输入，应该是系统不支持这种输入格式；所以必须每次运行都要包括build_orgdb，有点鸡肋。**Build_orgdb**可供需要建库个性分析使用

- **Image:**  
  Image: enrich-R--04; Enrich-R--06; go-figure

- **Results:**  
  

- **Build Database:**  
  ![build_orgdb_workflow](../png/build_orgdb_workflow.png)
  Ref
  > [模式植物构建orgDb数据库 | 以org.Slycompersicum.eg.db为例](https://mp.weixin.qq.com/s/b8OrDKJJGdXwF9B1C7l6zg)
  > [使用clusterProfiler对非模式植物进行注释](https://mp.weixin.qq.com/s/Mr3YLoc_-Y1WeLKJku1TzQ)
  > [富集分析|非模式物种GO/KEGG注释不会做？全网最详细eggNOG-mapper构建OrgDb包用于GO和KEGG富集分析](https://mp.weixin.qq.com/s/3sRdRuz6o5XuG11e2cX7Kw)
  > [生信干货 | AnnotationHub包-非模式物种OrgDB下载制作](https://mp.weixin.qq.com/s/auyTKJhfos0wi_yPsA7O0g)
  > [超详细非模式物种GO数据库全新代码构建](https://mp.weixin.qq.com/s/b23itzn5RNT8mJ1Ok8RzzA)

- **Using eggnog-mapper to Prepare Background Database Files** (related workflow: **Build_orgdb**)
  - [Recommended tutorial for building orgDb databases (in Chinese)](https://mp.weixin.qq.com/s/b8OrDKJJGdXwF9B1C7l6zg)
    1. Obtain the protein sequences for your species. Ideally, the protein names should match the gene names in your single-cell sequencing data
    2. Submit the protein sequences to [eggnog-mapper](http://eggnog-mapper.embl.de/): Select your protein sequence file; Enter your email address;Submit the job; Click the link in your email (`Click to manage your job`); Start the job and download the resulting `.xlsx` file.  
---

# Reference


---

# Coder info
- **Editor:** yangdong (yangdong@genomics.cn)
- **GitHub:** [ydgenomics](https://github.com/ydgenomics)
- **Prospect:** Focused on innovative, competitive, open-source projects.
- **Repository:** [Scripts/enrich_scRNAseq](https://github.com/ydgenomics/Scripts/tree/main/enrich_scRNAseq)