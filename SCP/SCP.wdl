version 1.0
workflow dataget_scRNAseqV1_2_1{
  input{
    Array[Array[File]] RawMatrix
    Array[Array[File]] FilterMatrix
    Array[Array[File]] SpliceMatrix
    Array[Array[File]] UnspliceMatrix
    Array[Array[String]] sample_value
    Array[String] biosample_value
    String species='zimia'
    Int mem_soupx=32
    Int mem_scrublet=16
    File? mitogenes_txt
    Int? mito_threshold
  }
  String soupx_env = "stereonote_hpc/yangdong_bca0ae6099a24097b3ddd7aabcb511b2_private:latest" #SoupX-R--02
  String dataget_env="stereonote_hpc/yangdong_1791d26bbea94753bd84206bccc75bab_private:latest" #scrublet-py--04
  String url_scdatacg="stereonote_hpc/yangdong_6c3a7cd28b5d4861ad87065a5644f7ca_private:latest"
  String s_env="stereonote_hpc/yangdong_fc53bb93d99a44d8b707bb9cf3f72bbc_private:latest"
  call wdl{
    input:
    env=s_env
  }
  Int jobn=length(biosample_value)
  scatter(index in range(jobn)){
    Array[File] RawMatrix0=RawMatrix[index]
    Array[File] FilterMatrix0=FilterMatrix[index]
    Array[File] SpliceMatrix0=SpliceMatrix[index]
    Array[File] UnspliceMatrix0=UnspliceMatrix[index]
    Array[String] Sample0=sample_value[index]
    call scrublet{
      input:
        Matrix=FilterMatrix0,
        SpliceMatrix=SpliceMatrix0,
        UnspliceMatrix=UnspliceMatrix0,
        sample=Sample0,
        species=biosample_value[index],
        group_key="sample",
        mitogenes_txt=mitogenes_txt,
        mito_threshold=mito_threshold,
        cpu=4,
        mem=mem_scrublet,
        env=dataget_env,
        mingenes=100,
        mincells=3,
    }
    Int jobn2=length(RawMatrix0)
    scatter(index2 in range(jobn2)){
      call soupx{
        input:
        RawMatrix=RawMatrix0[index2],
        FilterMatrix=FilterMatrix0[index2],
        script=wdl.r,
        sample=Sample0[index2],
        maxrho=0.2,
        cpu=8,
        mem=mem_soupx,
        env=soupx_env,
        minCG=100,
        tfidfMin=1,
      }
    }
    Array[File] SoupXMatrix=soupx.soupx_matrix
    Array[File] SoupXtxt=soupx.soupx_txt
    Array[File] SoupXpdf=soupx.soupx_pdf
    call scrublet as sscrublet{
      input:
        Matrix=SoupXMatrix,
        SpliceMatrix=SpliceMatrix0,
        UnspliceMatrix=UnspliceMatrix0,
        sample=Sample0,
        species=biosample_value[index]+"_soupx",
        group_key="sample",
        mitogenes_txt=mitogenes_txt,
        mito_threshold=mito_threshold,
        cpu=4,
        mem=mem_scrublet,
        env=dataget_env,
        mingenes=100,
        mincells=3,
        SoupXtxt=SoupXtxt,
        SoupXpdf=SoupXpdf,
    }
  }
  call wdl_result1{
    input:
    sx_files=sscrublet.endfile,
    unsx_files=scrublet.endfile,
    env=s_env
  }
  # ----------------------------------- enrich --------------------------------------
  # env
  String url_wdl_enrich_scRNAseq = "stereonote_hpc/yangdong_d40d88da7c6e47b79ff32104775dd643_private:latest"
  String url_enrich = "stereonote_hpc/yangdong_a59c2bdb0cc74697811f72af890fdfb9_private:latest" # enrich-R--04
  String url_gofigure = "stereonote_hpc/yangdong_62c01e5f5e724b32b2591c91a9a722e6_private:latest" # go-figure
  # local files
  File ko_json="/Files/yangdong/enrich/ko00001.json"
  File go_obo="/Files/yangdong/enrich/go.obo"
  File ic_tsv="/Files/yangdong/enrich/ic.tsv"
  File relations_full_tsv="/Files/yangdong/enrich/relations_full.tsv"
  String taxid="1111"
  Float minp=0.05
  call wdl_enrich_scRNAseq{
    input:
    url=url_wdl_enrich_scRNAseq,
  }
  call build_orgdb{
    input:
    emapper_xlsx=emapper_xlsx,
    ko_json=ko_json,
    go_obo=go_obo,
    genus=species,
    species=species,
    taxid=taxid,
    cpu=2,
    mem=16,
    url=url_enrich,
    script=wdl_enrich_scRNAseq.r1,
  }
  String url_marker_plot="stereonote_hpc/yangdong_53ba514f3a944bf5866dc7ac31358b9c_private:latest"
  String url_annotation="stereonote_hpc/yangdong_ad6cde50ecec4bfc87190a5b05cfd3c0_private:latest"
  Int jobn3=length(scrublet.h5ad)
  scatter(index3 in range(jobn3)){
    call enrich{
      input:
      DB=build_orgdb.db,
      gene_csv=scrublet.csv[index3],
      kegg_info_RData=build_orgdb.kegg_info_RData,
      genus=genus,
      species=biosample_value[index3],
      minp=minp,
      cpu=2,
      mem=8,
      url=url_enrich,
      script=wdl_enrich_scRNAseq.r2,
    }
    call gofigure{
      input:
      enrich_result=enrich.result,
      ic_tsv=ic_tsv,
      relations_full_tsv=relations_full_tsv,
      go_obo=go_obo,
      species=biosample_value[index3],
      max_label=15,
      cpu=2,
      mem=8,
      url=url_gofigure,
      script1=wdl_enrich_scRNAseq.r3,
      script2=wdl_enrich_scRNAseq.py1,
    }
  }
  Int jobn4=length(sscrublet.h5ad)
  scatter(index4 in range(jobn4)){
    call enrich as senrich{
      input:
      DB=build_orgdb.db,
      gene_csv=scrublet.csv[index4],
      kegg_info_RData=build_orgdb.kegg_info_RData,
      genus=genus,
      species=biosample_value[index4]+'_soupx',
      minp=minp,
      cpu=2,
      mem=8,
      url=url_enrich,
      script=wdl_enrich_scRNAseq.r2,
    }
    call gofigure as sgofigure{
      input:
      enrich_result=enrich.result,
      ic_tsv=ic_tsv,
      relations_full_tsv=relations_full_tsv,
      go_obo=go_obo,
      species=biosample_value[index4]+'_soupx',
      max_label=15,
      cpu=2,
      mem=8,
      url=url_gofigure,
      script1=wdl_enrich_scRNAseq.r3,
      script2=wdl_enrich_scRNAseq.py1,
    }
  }
  call wdl_result2{
    input:
    sx_enrich=senrich.result,
    unsx_enrich=enrich.result,
    sx_gofigure=sgofigure.result,
    unsx_gofigure=gofigure.result,
    env=s_env
  } 
  output{
    File result1=wdl_result1.result1
    File orgdb=build_orgdb.db
    File result2=wdl_result2.result2

  }
}

task wdl{
  input{
    String env
  }
  command <<<
    cp /Scripts/dataget_scRNAseq/run_SoupX.R run_SoupX.R
  >>>
  output{
    File r="run_SoupX.R"
  }
  runtime{
    docker_url: "~{env}"
    req_cpu: 1
    req_memory: "2Gi"
  }
}

task soupx{
  input{
    File RawMatrix
    File FilterMatrix
    File script
    String sample
    Float maxrho
    Int minCG=100
    Int tfidfMin=1
    String path="result"
    Int cpu
    Int mem
    String env
  }
  command <<<
    ###################### Input section ############################
    raw_path="~{RawMatrix}"
    filter_path="~{FilterMatrix}"
    sample_name="~{sample}"
    minCG=~{minCG}
    tfidfMin=~{tfidfMin}
    highestrho=~{maxrho}
    ################################################################
    /opt/conda/bin/Rscript ~{script} \
    --raw_path $raw_path --filter_path $filter_path --sample_name $sample_name \
    --minCG $minCG --tfidfMin $tfidfMin --highestrho $highestrho
  >>>
  output{
    File soupx_txt="${sample}_soupx_rho.txt"
    File soupx_pdf="${sample}_rho.pdf"
    File soupx_matrix="${sample}"
  }
  runtime{
    docker_url: "~{env}"
    req_cpu: cpu
    req_memory: "~{mem}Gi"
  }
}

task scrublet{
  input{
    Array[File] Matrix #FilterMatrix or SoupxMatrix
    Array[File] SpliceMatrix
    Array[File] UnspliceMatrix
    Array[String] sample
    String species
    String group_key
    Int mingenes
    Int mincells
    File? mitogenes_txt
    Int? mito_threshold
    Int cpu
    Int mem
    String env
    Array[File]? SoupXtxt
    Array[File]? SoupXpdf
  }
  String outfile=species+"_dataget"
  command <<<
    mkdir ~{outfile}
    cd ~{outfile}
    for c in ~{sep="," Matrix}; do
        echo $c >> Matrix.txt
    done
    for c in ~{sep="," SpliceMatrix}; do
        echo $c >> SpliceMatrix.txt
    done
    for c in ~{sep="," UnspliceMatrix}; do
        echo $c >> UnspliceMatrix.txt
    done
    for d in ~{sep="," sample}; do
        echo $d >> samples.txt
    done
    for c in $(echo ~{sep="," SoupXtxt} | tr ',' ' '); do
        echo "$c"
        cp $c .
    done
    for c in $(echo ~{sep="," SoupXpdf} | tr ',' ' '); do
        echo "$c"
        cp $c .
    done
    samples_txt_path="samples.txt"
    /opt/conda/bin/python << CODE
    species = "~{species}"
    group_key = "~{group_key}"
    matrix_txt = "Matrix.txt"
    splice_txt = "SpliceMatrix.txt"
    unsplice_txt = "UnspliceMatrix.txt"
    sample_txt = "samples.txt"
    input_mingenes = ~{mingenes}
    input_mincells = ~{mincells}
    mito_genes = "~{mitogenes_txt}"
    if mito_genes == "":
        mito_genes = "None_mito_genes.csv"
    
    #mito_threshold = ~{mito_threshold}
    try:
        mito_threshold = float("~{mito_threshold}") 
    except ValueError:
        mito_threshold = 0.05
    
    # Date: 20250607
    # Attention: how to rationally get a multi-matrix anndata including FilterMatrix, SpliceMatrix and UnspliceMatrix.
    # Marker_csv: gene, cluster, p_val_adj, avg_log2FC

    import numpy as np
    import pandas as pd
    import scanpy as sc
    import anndata as ad
    import seaborn as sns
    from matplotlib.pyplot import savefig
    from pathlib import Path
    import shutil
    import gzip
    import os
    import sys
    import scrublet
    import leidenalg
    import argparse
    # import logging

    # Get command line arguments
    parser = argparse.ArgumentParser(description="Estimate double cells using Scrublet and process multi-matrix AnnData.")
    parser.add_argument('--species', type=str, default='zimia', help='Species name')
    parser.add_argument('--group_key', type=str, default='sample', help='Group key for batch')
    parser.add_argument('--matrix_txt', type=str, default="Matrix.txt", help='Path to matrix file list')
    parser.add_argument('--splice_txt', type=str, default="SpliceMatrix.txt", help='Path to splice file list')
    parser.add_argument('--unsplice_txt', type=str, default="UnspliceMatrix.txt", help='Path to unsplice file list')
    parser.add_argument('--sample_txt', type=str, default="samples.txt", help='Path to sample names file')
    parser.add_argument('--input_mingenes', type=int, default=100, help='Minimum number of genes per cell')
    parser.add_argument('--input_mincells', type=int, default=3, help='Minimum number of cells per gene')
    parser.add_argument('--mito_genes', type=str, default="None_mito_genes.csv", help='CSV file with mitochondrial genes')
    parser.add_argument('--mito_threshold', type=float, default=0.05, help='Mitochondrial gene threshold')

    args = parser.parse_args()
    
    def copy_and_process(matrixfile, featuresfile, barcodesfile, target_folder):
        """
        Copy and decompress matrix, features, and barcodes files to the target folder.
        """
        original_dir = os.getcwd()
        os.chdir(target_folder)
        shutil.copy(matrixfile, "matrix.mtx.gz")
        shutil.copy(featuresfile, "features.tsv.gz")
        shutil.copy(barcodesfile, "barcodes.tsv.gz")
        with gzip.open('matrix.mtx.gz', 'rb') as g_file1, open("matrix.mtx", "wb") as f_out:
            f_out.write(g_file1.read())
        with gzip.open('features.tsv.gz', 'rb') as g_file2, open("features.tsv", "wb") as f_out:
            f_out.write(g_file2.read())
        with gzip.open('barcodes.tsv.gz', 'rb') as g_file3, open("barcodes.tsv", "wb") as f_out:
            f_out.write(g_file3.read())
        with open('features.tsv', 'r') as f_in, open('genes.tsv', 'w') as f_out:
            for line in f_in:
                f_out.write(line.strip() + '\t' + line.strip() + '\n')
        os.chdir(original_dir)

    def complete_genes(adata, all_genes, gene_symbols_col='gene_symbols'):
        """
        Complete missing genes in the AnnData object and set their values to 0.

        Args:
            adata (AnnData): AnnData object to be completed.
            all_genes (set): Complete set of genes.
            gene_symbols_col (str): Column name for gene symbols, default is 'gene_symbols'.

        Returns:
            AnnData: AnnData object with completed genes.
        """
        current_genes = set(adata.var_names)
        missing_genes = all_genes - current_genes

        if len(missing_genes) > 0:
            print(f"Completing missing genes: {len(missing_genes)}")
            missing_genes_df = pd.DataFrame(
                0, index=adata.obs_names, columns=list(missing_genes)
            )
            missing_genes_adata = ad.AnnData(
                X=missing_genes_df.values,
                obs=adata.obs,
                var=pd.DataFrame(index=list(missing_genes))
            )
            missing_genes_adata.var[gene_symbols_col] = missing_genes_adata.var.index
            adata = ad.concat([adata, missing_genes_adata], axis=1)
            adata = adata[:, list(all_genes)]
        else:
            print("No need to complete, all genes are present in adata.")
        return adata

    def complete_cells(adata, all_cells):
        """
        Complete missing cells in the AnnData object and set their values to 0.

        Args:
            adata (AnnData): AnnData object to be completed.
            all_cells (set): Complete set of cells.

        Returns:
            AnnData: AnnData object with completed cells.
        """
        current_cells = set(adata.obs_names)
        missing_cells = all_cells - current_cells

        if len(missing_cells) > 0:
            print(f"Completing missing cells: {len(missing_cells)}")
            missing_cells_df = pd.DataFrame(
                0, index=list(missing_cells), columns=adata.var_names
            )
            missing_cells_adata = ad.AnnData(
                X=missing_cells_df.values,
                obs=pd.DataFrame(index=list(missing_cells)),
                var=adata.var
            )
            adata = ad.concat([adata, missing_cells_adata], axis=0)
            adata = adata[list(all_cells), :]
        else:
            print("No need to complete, all cells are present in adata.")
        return adata

    def run_concat_plot(species, input_mingenes, input_mincells, group_key, sample_names, trans_matrix_list, trans_splice_list, trans_unsplice_list, mito_genes, mito_threshold):
        """
        Concatenate, QC, filter, and plot AnnData objects for all samples.
        """
        adatas = {}
        for i in range(len(sample_names)): 
            key = sample_names[i]
            adata_filter = sc.read_10x_mtx(trans_matrix_list[i], var_names='gene_ids')
            adata_splice = sc.read_10x_mtx(trans_splice_list[i], var_names='gene_ids')
            adata_unsplice = sc.read_10x_mtx(trans_unsplice_list[i], var_names='gene_ids')
            # Get gene sets for each dataset
            genes_filter = set(adata_filter.var_names)
            genes_splice = set(adata_splice.var_names)
            genes_unsplice = set(adata_unsplice.var_names)
            all_genes = genes_filter.union(genes_splice).union(genes_unsplice)
            print(f"sample: {key}, genes in matrix/splice/unsplice/union: {len(genes_filter)}/{len(genes_splice)}/{len(genes_unsplice)}/{len(all_genes)}")
            adata_filter = complete_genes(adata_filter, all_genes)
            adata_splice = complete_genes(adata_splice, all_genes)
            adata_unsplice = complete_genes(adata_unsplice, all_genes)
            # Get cell sets for each dataset
            cells_filter = set(adata_filter.obs_names)
            cells_splice = set(adata_splice.obs_names)
            cells_unsplice = set(adata_unsplice.obs_names)
            all_cells = cells_filter.union(cells_splice).union(cells_unsplice)
            print(f"sample: {key}, cells in matrix/splice/unsplice/union: {len(cells_filter)}/{len(cells_splice)}/{len(cells_unsplice)}/{len(all_cells)}")
            adata_filter = complete_cells(adata_filter, all_cells)
            adata_splice = complete_cells(adata_splice, all_cells)
            adata_unsplice = complete_cells(adata_unsplice, all_cells)
            adata = adata_filter.copy()
            adata.layers['splice'] = adata_splice.X
            adata.layers['unsplice'] = adata_unsplice.X
            # Rename cells to include sample key
            adata.obs_names = [f"{cell_name}_{key}" for cell_name in adata.obs_names]
            print(adata.obs_names[:10])
            adatas[key] = adata
        adata = ad.concat(adatas, label=group_key, join="outer")
        adata.obs['biosample'] = species
        print(adata.obs[group_key].value_counts())

        # Set parameters for figures
        sc.settings.verbosity = 3
        sc.logging.print_versions()
        sc.settings.set_figure_params(dpi=80, facecolor='white')

        # Check mitochondrial genes and filter
        if os.path.exists(mito_genes):
            mt_genes = pd.read_csv(mito_genes, header=None, names=["gene_name"])
            mt_genes_list = mt_genes["gene_name"].tolist()
            print(mt_genes_list[:10])
            adata.var["mt"] = adata.var_names.isin(mt_genes)
            print("calculate mt genes")
            sc.pp.calculate_qc_metrics(adata,qc_vars=["mt"],inplace=True,log1p=True)
            sc.pl.violin(adata,["n_genes_by_counts", "total_counts", "pct_counts_mt"],jitter=0.4,multi_panel=True,save="_mitogene.pdf")
            sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", save="_mitogenes.pdf")
            adata = adata[adata.obs.pct_counts_mt < mito_threshold].copy()
            sc.pl.violin(adata,["n_genes_by_counts", "total_counts", "pct_counts_mt"],jitter=0.4,multi_panel=True,save="_mitogene_filtered.pdf")
            sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", save="_mitogenes_filtered.pdf")
        else:
            print("mitochondrial list not exist")
            sc.pp.calculate_qc_metrics(adata, inplace=True, log1p=True)
        sns.jointplot(data=adata.obs, x="log1p_total_counts", y="log1p_n_genes_by_counts", kind="hex")
        savefig("qc.pdf")

        # Pre-process, QC, and Scrublet
        sc.pp.filter_cells(adata, min_genes=input_mingenes)
        sc.pp.filter_genes(adata, min_cells=input_mincells)
        sc.external.pp.scrublet(adata, batch_key=group_key)

        adata.layers["counts"] = adata.X.copy()

        # Visualization
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key=group_key)
        sc.tl.pca(adata)
        # Check if group_key exists in obs
        if group_key not in adata.obs:
            raise ValueError(f"Group key '{group_key}' not found in adata.obs")
        features = [group_key, group_key]
        if 'pct_counts_mt' in adata.obs:
            features.extend(['pct_counts_mt', 'pct_counts_mt'])
        features.extend(['doublet_score', 'doublet_score'])
        dimensions = [(0, 1), (2, 3)] * (len(features) // 2)
        save_filename = '_potentially_undesired_features'
        if 'pct_counts_mt' in adata.obs:
            save_filename += '_with_mt'
        save_filename += '.pdf'
        sc.pl.pca(adata, color=features, dimensions=dimensions, ncols=2, size=2, save=save_filename)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        sc.pl.umap(adata, color=group_key, size=2, save="_batch.pdf")
        sc.tl.leiden(adata, resolution=1)
        adata.obs['predicted_doublet'] = adata.obs['predicted_doublet'].astype('category')
        sc.pl.umap(adata, color=["leiden", "log1p_n_genes_by_counts", "predicted_doublet", "doublet_score"], ncols=2, save="_quality.pdf")
        for res in [0.02, 0.2, 0.5, 0.8, 1.0, 1.3, 1.6, 2.0]:
            sc.tl.leiden(adata, key_added=f"leiden_res_{res:4.2f}", resolution=res)
        sc.pl.umap(adata, color=["leiden_res_0.02", "leiden_res_0.20", "leiden_res_0.50", "leiden_res_0.80", "leiden_res_1.00", "leiden_res_1.30", "leiden_res_1.60", "leiden_res_2.00"], legend_loc="on data", save="_leiden_clus.pdf")
        # Marker
        output_dir = "marker_csv"
        os.makedirs(output_dir)
        resolutions = ["leiden_res_0.50", "leiden_res_0.80", "leiden_res_1.00"]
        for res in resolutions:
            sc.tl.rank_genes_groups(adata, groupby=res, method="wilcoxon")
            sc.pl.rank_genes_groups_dotplot(adata, groupby=res, standard_scale="var", n_genes=5, save=f"{res}_marker.pdf")
            marker = sc.get.rank_genes_groups_df(adata, group=None)
            marker['gene'] = marker['names']
            marker['cluster'] = marker['group']
            marker['p_val_adj'] = marker['pvals_adj']
            marker['avg_log2FC'] = marker['logfoldchanges']
            marker.to_csv(f"{output_dir}/{res}.markers.csv", index=False)
        # Summary
        with open('summary.txt', 'w') as f:
            f.write(species + ' data summary' + '\n')
            f.write('Total cells: ' + str(adata.n_obs) + '\n')
            f.write('Total genes: ' + str(adata.n_vars) + '\n')
            f.write('Average genes per cell: ' + str(adata.obs['n_genes'].mean()) + '\n')
            f.write('Median genes per cell: ' + str(adata.obs['n_genes'].median()) + '\n')
            f.write('Average counts per cell: ' + str(adata.obs['total_counts'].mean()) + '\n')
            f.write('Median counts per cell: ' + str(adata.obs['total_counts'].median()) + '\n')
            # Write top 10 cell and gene names
            f.write('\nTop 10 cells:\n' + ','.join(adata.obs_names[:10]) + '\n')
            f.write('\nTop 10 genes:\n' + ','.join(adata.var_names[:10]) + '\n')
        adata.X = adata.layers["counts"] # Save the raw counts in the X attribute
        adata.write_h5ad(filename=species + '.h5ad', compression="gzip")

    # Main function to run the scrublet analysis
    def run_scrublet(species, sample_txt, matrix_txt, splice_txt, unsplice_txt, input_mingenes=100, input_mincells=3, group_key="sample", mito_genes="None_mito_genes.csv", mito_threshold=0.05):
        """
        Main function to run Scrublet and process multi-matrix AnnData.
        """
        # Load the data: from text transform to array
        with open(matrix_txt, 'r') as file:
            matrix_files = file.read().strip().split(',')
        with open(splice_txt, 'r') as file:
            splice_files = file.read().strip().split(',')
        with open(unsplice_txt, 'r') as file:
            unsplice_files = file.read().strip().split(',')
        with open(sample_txt, 'r') as filen:
            sample_names = filen.read().strip().split(',')

        # Preprocess the loaded data
        trans_matrix_list = []
        trans_splice_list = []
        trans_unsplice_list = []

        process_types = [
            ("filter", matrix_files),
            ("splice", splice_files),
            ("unsplice", unsplice_files)
        ]
        if len(matrix_files) > 0:
            for i in range(len(sample_names)):
                sample = sample_names[i]
                for process_name, file_list in process_types:
                    directory_path = Path(f"./{sample}/{process_name}")
                    directory_path.mkdir(parents=True, exist_ok=True)
                    folder_path = os.path.abspath(directory_path)
                    if process_name == "filter":
                        trans_matrix_list.append(folder_path)
                        matrixfile = file_list[i] + '/matrix.mtx.gz'
                        featuresfile = file_list[i] + '/features.tsv.gz'
                        barcodesfile = file_list[i] + '/barcodes.tsv.gz'
                    elif process_name == "splice":
                        trans_splice_list.append(folder_path)
                        matrixfile = file_list[i] + '/matrix.mtx.gz'
                        featuresfile = file_list[i] + '/features.tsv.gz'
                        barcodesfile = file_list[i] + '/barcodes.tsv.gz'
                    elif process_name == "unsplice":
                        trans_unsplice_list.append(folder_path)
                        matrixfile = file_list[i] + '/unspliced.mtx.gz'
                        featuresfile = file_list[i] + '/features.tsv.gz'
                        barcodesfile = file_list[i] + '/barcodes.tsv.gz'
                    copy_and_process(matrixfile, featuresfile, barcodesfile, folder_path)
            print(trans_matrix_list); print(trans_splice_list); print(trans_unsplice_list); print(sample_names)
            run_concat_plot(species, input_mingenes, input_mincells, group_key, sample_names, trans_matrix_list, trans_splice_list, trans_unsplice_list, mito_genes, mito_threshold)
        else:
            print("No samples to process")
            with open('summary.txt', 'w') as f:
                f.write(species + ' data summary' + '\n')
                f.write('No samples to process' + '\n')
    
    run_scrublet(species, sample_txt, matrix_txt, splice_txt, unsplice_txt, input_mingenes, input_mincells, group_key, mito_genes, mito_threshold)
    CODE
    file_paths=$(tr ',' '\n' < $samples_txt_path)
    echo $file_paths
    for path in ${file_paths[@]}; do
        rm -rf $path
    done
    rm "Matrix.txt"
    rm "SpliceMatrix.txt"
    rm "UnspliceMatrix.txt"
    rm "samples.txt"
  >>>
  output{
    File endfile="~{outfile}"
    File h5ad="~{outfile}/~{species}.h5ad"
    File csv="~{outfile}/marker_csv/leiden_res_0.50.markers.csv"
  }
  runtime{
    docker_url: "~{env}"
    req_cpu: cpu
    req_memory: "~{mem}Gi"
  }
}
task wdl_result1{
  input{
    Array[File] sx_files
    Array[File] unsx_files
    String env
  }
  command <<<
    mkdir 01_dataget
    cd 01_dataget
    for c in $(echo ~{sep="," sx_files} | tr ',' ' '); do
        echo "$c"
        input_folder="$c"
        folder_name=$(basename "$input_folder")
        tar -czvf "$folder_name".tar.gz -C "$(dirname "$input_folder")" "$folder_name"
        tar -zxvf "$folder_name".tar.gz
    done
    for c in $(echo ~{sep="," unsx_files} | tr ',' ' '); do
        echo "$c"
        input_folder="$c"
        folder_name=$(basename "$input_folder")
        tar -czvf "$folder_name".tar.gz -C "$(dirname "$input_folder")" "$folder_name"
        tar -zxvf "$folder_name".tar.gz
    done
  >>>
  output{
    File result1="01_dataget"
  }
  runtime{
    docker_url: "~{env}"
    req_cpu: 1
    req_memory: "4Gi"
  }
}

task wdl_enrich_scRNAseq{
  input{
    String url
  }
  command <<<
    cp /Scripts/enrich_scRNAseq/build_orgdb.R build_orgdb.R
    cp /Scripts/enrich_scRNAseq/enrich.R enrich.R
    cp /Scripts/enrich_scRNAseq/deal_enrich_txt.R deal_enrich_txt.R
    cp /Scripts/enrich_scRNAseq/gofigure.py gofigure.py
  >>>
  output{
    File r1="build_orgdb.R"
    File r2="enrich.R"
    File r3="deal_enrich_txt.R"
    File py1="gofigure.py"
  }
  runtime{
    docker_url: "~{url}"
    req_cpu: 1
    req_memory: "2Gi"
  }
}

task build_orgdb{
  input {
    File emapper_xlsx
    File ko_json
    File go_obo
    String genus
    String species
    String taxid
    Int cpu
    Int mem
    String url
    File script
  }
  command <<<
    mkdir build_orgdb
    cd build_orgdb
    
    emapper_xlsx="~{emapper_xlsx}"
    ko_json="~{ko_json}"
    go_obo="~{go_obo}"
    taxid="~{taxid}"
    genus="~{genus}"
    species="~{species}"
    
    /opt/conda/bin/Rscript ~{script} \
    --emapper_xlsx $emapper_xlsx --ko_json $ko_json --go_obo $go_obo \
    --taxid $taxid --genus $genus --species $species
  >>>
  runtime {
    docker_url: "~{url}"
    req_cpu: cpu
    req_memory: "~{mem}Gi"
  }
  output {
    File db = "build_orgdb"
    File kegg_info_RData = "build_orgdb/kegg_info.RData"
  }
}

task enrich{
  input {
    File DB
    File gene_csv
    File kegg_info_RData
    String genus
    String species
    Float minp
    Int cpu
    Int mem
    String url
    File script
  }
  command <<<
    gene_csv="~{gene_csv}"
    kegg_info_RData="~{kegg_info_RData}"
    db="~{DB}"
    minp=~{minp}
    genus="~{genus}"
    species="~{species}"
    
    /opt/conda/bin/Rscript ~{script} \
    --gene_csv $gene_csv --kegg_info_RData $kegg_info_RData --db $db \
    --minp $minp --genus $genus --species $species
  >>>
  runtime {
    docker_url: "~{url}"
    req_cpu: cpu
    req_memory: "~{mem}Gi"
  }
  output {
    File result="${species}_enrich"
  }
}

task gofigure{
  input {
    File enrich_result
    File ic_tsv
    File relations_full_tsv
    File go_obo
    String species
    Int max_label
    Int cpu
    Int mem
    String url
    File script1
    File script2
  }
  command <<<
    /opt/software/R/bin/Rscript ~{script1} ~{enrich_result}

    mkdir "${species}_gofigure_result"
    cd "${species}_gofigure_result"
    cp ~{script2} .
    
    result_name_file="../result_name.txt"
    output_file="../output_standard_gofigure_input.txt"
    max_label=~{max_label}

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
    # run go-figure
    mkdir "data"
    cp ~{ic_tsv} "data"
    cp ~{relations_full_tsv} "data"
    cp ~{go_obo} "data"
    # 循环处理每个名称和对应的输出文件
    i=0
    while [ $i -lt $len_names ]; do
      name=${names[$i]}
      output=${outputs[$i]}

      echo "Processing $name with output file $output"
      mkdir "$name"
      tsv_path="$output"
      max_label=15

      mkdir "$name"
      /software/miniconda/envs/go-figure/bin/python gofigure.py \
        -i "$tsv_path" -j standard -m "$max_label" -o "$name"
      # 自增索引
      i=$((i + 1))
    done
    rm -r "data"
    rm gofigure.py
  >>>
  runtime {
    docker_url: "~{url}"
    req_cpu: cpu
    req_memory: "~{mem}Gi"
  }
  output {
    File result="${species}_gofigure_result
  }
}

task wdl_result2{
  input{
    Array[File] sx_enrich
    Array[File] unsx_enrich
    Array[File] sx_gofigure
    Array[File] unsx_gofigure
    String env
  }
  command <<<
    mkdir 02_enrich
    cd 02_enrich
    for c in $(echo ~{sep="," sx_enrich} | tr ',' ' '); do
        echo "$c"
        input_folder="$c"
        folder_name=$(basename "$input_folder")
        tar -czvf "$folder_name".tar.gz -C "$(dirname "$input_folder")" "$folder_name"
        tar -zxvf "$folder_name".tar.gz
    done
    for c in $(echo ~{sep="," unsx_enrich} | tr ',' ' '); do
        echo "$c"
        input_folder="$c"
        folder_name=$(basename "$input_folder")
        tar -czvf "$folder_name".tar.gz -C "$(dirname "$input_folder")" "$folder_name"
        tar -zxvf "$folder_name".tar.gz
    done
    for c in $(echo ~{sep="," sx_gofigure} | tr ',' ' '); do
        echo "$c"
        input_folder="$c"
        folder_name=$(basename "$input_folder")
        tar -czvf "$folder_name".tar.gz -C "$(dirname "$input_folder")" "$folder_name"
        tar -zxvf "$folder_name".tar.gz
    done
    for c in $(echo ~{sep="," unsx_gofigure} | tr ',' ' '); do
        echo "$c"
        input_folder="$c"
        folder_name=$(basename "$input_folder")
        tar -czvf "$folder_name".tar.gz -C "$(dirname "$input_folder")" "$folder_name"
        tar -zxvf "$folder_name".tar.gz
    done
  >>>
  output{
    File result2="02_enrich"
  }
  runtime{
    docker_url: "~{env}"
    req_cpu: 1
    req_memory: "4Gi"
  }
}

task check_rds{
  input {
    File rds
    String save_rds
    Int mem
    String url
  }
  command {
    /opt/conda/bin/R --vanilla --slave <<EOF
    library(Seurat)
    seu <- readRDS("~{rds}")
    print("####################################")
    print(seu)
    print("####################################")
    print(colnames(seu@meta.data))
    print("####################################")
    print(head(seu@meta.data))
    
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
    seu <- ScaleData(seu)
    
    saveRDS(seu, file="~{save_rds}")
    EOF
  }
  runtime {
    docker_url: "~{url}"
    req_cpu: 2
    req_memory: "~{mem}Gi"
  }
  output {
    File checked_rds = "~{save_rds}"
  }
}

task marker_plot{
  input {
    File rds
    File markers_csv
    File cluster_color_csv
    String tissue
    String cluster_key
    String reduction_key
    String species
    Int mem
    String url
  }
  command {
    mkdir "${species}_marker_plot"
    cd "${species}_marker_plot"
    
    input_rds="~{rds}"
    markers_csv="~{markers_csv}"
    cluster_color_csv="~{cluster_color_csv}"
    cell_type="~{tissue}"
    cluster_key="~{cluster_key}"
    reduction_key="~{reduction_key}"

    /opt/conda/bin/Rscript /script/visual_cg.R \
    --input_rds $input_rds --markers_csv $markers_csv \
    --cluster_color_csv $cluster_color_csv --cell_type $cell_type \
    --cluster_key $cluster_key --reduction_key $reduction_key
  }
  runtime {
    docker_url: "~{url}"
    req_cpu: 2
    req_memory: "~{mem}Gi"
  }
  output {
    File plot = "${species}_marker_plot"
  }
}

task sctype{
  input {
    File input_query_rds
    File input_marker_csv
    String tissue
    String cluster_key
    String umap_name
    String species
    Int n_circle
    Int mem
    String url
  }
  command {
    mkdir "${species}_anno_sctype"
    cd "${species}_anno_sctype"
    
    input_query_rds="~{input_query_rds}"
    input_marker_csv="~{input_marker_csv}"
    tissue="~{tissue}"
    cluster_key="~{cluster_key}"
    umap_name="~{umap_name}"
    n_circle=~{n_circle}

    /software/miniconda/envs/Seurat/bin/Rscript /script/anno/anno_sctype.R \
    --input_query_rds $input_query_rds --input_marker_csv $input_marker_csv \
    --tissue $tissue --n_circle $n_circle \
    --cluster_key $cluster_key --umap_name $umap_name
  }
  runtime {
    docker_url: "~{url}"
    req_cpu: 2
    req_memory: "~{mem}Gi"
  }
  output {
    File result=""${species}_anno_sctype""
  }
}