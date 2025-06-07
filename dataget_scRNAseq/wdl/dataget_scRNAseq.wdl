version 1.0
workflow dataget_scRNAseqV1_2_0{
  input{
    Array[File] RawMatrix
    Array[File] FilterMatrix
    Array[File] SpliceMatrix
    Array[File] UnspliceMatrix
    Array[File] Sample
    String Species
    Float Maxrho=0.2
    String Group_key="sample"
    Int mem_soupx=20
    Int mem_scrublet=10
    Int mem_scdatacg=8
    File? mitogenes_txt
    Int? mito_threshold
  }
  String soupx_env = "stereonote_hpc/yangdong_929e6c6e283247d28d4dd4d56d9c868d_private:latest" #SoupX-R--03
  String dataget_env="stereonote_hpc/yangdong_1791d26bbea94753bd84206bccc75bab_private:latest" #scrublet-py--04
  String url_scdatacg="stereonote_hpc/yangdong_6c3a7cd28b5d4861ad87065a5644f7ca_private:latest"
  call scrublet{
    input:
      Matrix=FilterMatrix,
      SpliceMatrix=SpliceMatrix,
      UnspliceMatrix=UnspliceMatrix,
      sample=Sample,
      species=Species,
      group_key=Group_key,
      mitogenes_txt=mitogenes_txt,
      mito_threshold=mito_threshold,
      cpu=4,
      mem=mem_scrublet,
      env=dataget_env,
      mingenes=100,
      mincells=3,
  }
  Int jobn=length(RawMatrix)
  scatter(index in range(jobn)){
    call soupx{
      input:
      RawMatrix=RawMatrix[index],
      FilterMatrix=FilterMatrix[index],
      sample=Sample[index],
      maxrho=Maxrho,
      cpu=10,
      mem=mem_soupx,
      env=soupx_env,
      minCG=100,
      tfidfMin=1,
    }
  }
  Array[File] SoupXMatrix=soupx.soupx_matrix
  call scrublet as sscrublet{
    input:
      Matrix=SoupXMatrix,
      SpliceMatrix=SpliceMatrix,
      UnspliceMatrix=UnspliceMatrix,
      sample=Sample,
      species=Species+"_soupx",
      group_key=Group_key,
      mitogenes_txt=mitogenes_txt,
      mito_threshold=mito_threshold,
      cpu=4,
      mem=mem_scrublet,
      env=dataget_env,
      mingenes=100,
      mincells=3,
  }
  Array[File] all_files=[scrublet.h5ad, sscrublet.h5ad]
  Int jobn2=length(all_files)
  scatter(index2 in range(jobn2)){
    call scdatacg{
      input:
      rds_h5ad=all_files[index2],
      layers='all',
      cpu=2,
      mem=mem_scdatacg,
      url=url_scdatacg,
    }
  }
  output{
    File outdir=scrublet.endfile
    File soutdir=sscrublet.endfile
    Array[File] rds=scdatacg.rds
    Array[File] rho_txt=soupx.soupx_txt
    Array[File] rho_pdf=soupx.soupx_pdf
  }
}

task soupx{
  input{
    File RawMatrix
    File FilterMatrix
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
    /opt/conda/bin/Rscript /script/dataget_scRNAseq/V1.2.0/run_SoupX.R \
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
  }
  runtime{
    docker_url: "~{env}"
    req_cpu: cpu
    req_memory: "~{mem}Gi"
  }
}

task scdatacg{
  input {
    File rds_h5ad
    String layers
    Int cpu
    Int mem
    String url
  }
  command <<<
    #!/bin/bash
    input_file=~{rds_h5ad}
    ext="${input_file##*.}"
    echo "input file extension is: $ext"
    if [ "$ext" == "rds" ]; then
        echo "Converting rds to h5ad..."
        /software/conda/Anaconda/bin/Rscript /script/convert_rdsAh5ad2.R --input_file $input_file --layers ~{layers}
        python /script/deal_layers_ydgenomics.py --input_path $input_file --sctype $ext
        echo "Copying rds file..."
        cp "$input_file" ./
    elif [ "$ext" == "h5ad" ]; then
        echo "Converting h5ad to rds..."
        python /script/deal_layers_ydgenomics.py --input_path $input_file --sctype $ext
        /software/conda/Anaconda/bin/Rscript /script/convert_rdsAh5ad2.R --input_file $input_file --layers ~{layers}
    else
        echo "Error: Unsupported file extension '$ext'. Only 'rds' and 'h5ad' are supported."
    fi
  >>>
  runtime {
    docker_url: "~{url}"
    req_cpu: cpu
    req_memory: "~{mem}Gi" 
  }
  output {
    File rds = glob("*.rds")[0]
  }
}