# /software/conda/Anaconda/bin/python 
# scib_metrics_benchmarking_test.py
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import webbrowser
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
import click

@click.command()
@click.argument("unintegrated_h5ad", type=click.Path(exists=True))
@click.argument("integrated_file", type=click.Path(exists=True))
@click.argument("methods_file", type=click.Path(exists=True))
@click.argument("pcas_file", type=click.Path(exists=True))
@click.argument("deals_file", type=click.Path(exists=True))
@click.argument("tests_file", type=click.Path(exists=True))
@click.option('--batch_key', type=str, default=None, help="Batch key")
@click.option('--label_key', type=str, default=None, help="Storying the information of biological cell name")
@click.option('--cluster_key', type=str, default=None, help="Cluster key in data to use as labels to marker different clusters")
@click.option('--n_jobs', type=int, default=None, help="Number of jobs to use for parallelization of neighbor search")
@click.argument("out_rawpdf", type=click.Path(exists=False), default=None)
@click.argument("out_benchpdf", type=click.Path(exists=False), default="benchmark_results.pdf")
@click.argument("out_benchcsv", type=click.Path(exists=False), default="benchmark_results.csv")
def main(unintegrated_h5ad, integrated_file, methods_file, pcas_file, deals_file, tests_file, batch_key, label_key, cluster_key, n_jobs, out_rawpdf, out_benchpdf, out_benchcsv):
    # Read files and split by comma
    with open(integrated_file, 'r') as file:
        files = file.read().strip().split(',')
    with open(methods_file, 'r') as file:
        methods = file.read().strip().split(',')
    with open(pcas_file, 'r') as file:
        pcas = file.read().strip().split(',')
    with open(deals_file, 'r') as file:
        deals = file.read().strip().split(',')
    with open(tests_file, 'r') as file:
        tests = file.read().strip().split(',')
    print(tests)
    
    # Process unintegrated data
    orig_ad = sc.read_h5ad(unintegrated_h5ad)
    sc.set_figure_params(dpi_save=200, frameon=False, figsize=(10, 5))
    sc.pp.normalize_total(orig_ad, target_sum=1e4)
    sc.pp.log1p(orig_ad)
    sc.pp.highly_variable_genes(orig_ad, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.scale(orig_ad, max_value=10)
    sc.tl.pca(orig_ad, svd_solver="arpack")
    sc.pp.neighbors(orig_ad, n_neighbors=20, n_pcs=20)
    sc.tl.leiden(orig_ad, resolution=0.5, key_added=cluster_key, flavor='igraph', n_iterations=2, directed=False)
    sc.tl.umap(orig_ad, min_dist=0.3)
    sc.pl.umap(orig_ad, color=[batch_key, label_key, cluster_key])
    plt.savefig(out_rawpdf, dpi=300, bbox_inches='tight')
    orig_ad.obsm["Unintegrated"] = orig_ad.obsm["X_pca"]

    # Process integrated data and merge with unintegrated data
    for i in range(len(pcas)):
        adata = sc.read_h5ad(files[i])
        if deals[i] == 'N':
            adata.obsm[methods[i]] = adata.obsm[pcas[i]]
            orig_ad.obsm[methods[i]] = adata.obsm[methods[i]]
        else:
            # deals[i] == 'Y', so we need to process the data which mostly have been processed by R methods
            if '_index' in adata.obs:
                del adata.obs['_index']
            if '_index' in adata.var:
                del adata.var['_index']
            sc.tl.pca(adata, svd_solver="arpack")
            adata.obsm[pcas[i]] = adata.obsm["X_pca"]
            adata.obsm[methods[i]] = adata.obsm[pcas[i]]
            sc.pp.neighbors(adata, n_neighbors=20, n_pcs=20, use_rep=pcas[i])
            orig_ad.obsm[methods[i]] = adata.obsm[methods[i]]

    #adata.write("merged_obsm_ad.h5ad", compression="gzip")

    methods.append('Unintegrated')

    import time
    start = time.time()
    def str_to_bool(value):
        return value.lower() in ("true", "yes", "1", "on")
    biocons = BioConservation(isolated_labels=str_to_bool(tests[0]), nmi_ari_cluster_labels_leiden=str_to_bool(tests[1]), nmi_ari_cluster_labels_kmeans=str_to_bool(tests[2]), silhouette_label=str_to_bool(tests[3]), clisi_knn=str_to_bool(tests[4]))
    bacorrec = BatchCorrection(silhouette_batch=str_to_bool(tests[5]), ilisi_knn=str_to_bool(tests[6]), kbet_per_label=str_to_bool(tests[7]), graph_connectivity=str_to_bool(tests[8]), pcr_comparison=str_to_bool(tests[9]))
    bm = Benchmarker(
        adata=orig_ad,
        batch_key=batch_key,
        label_key=label_key,
        embedding_obsm_keys=methods,
        pre_integrated_embedding_obsm_key="X_pca",
        bio_conservation_metrics=biocons,
        batch_correction_metrics=bacorrec,
        n_jobs=n_jobs,
    )
    bm.benchmark()
    end = time.time()
    print(f"Time: {int((end - start) / 60)} min {int((end - start) % 60)} sec")

    bm.plot_results_table()
    bm.plot_results_table(min_max_scale=False)
    plt.savefig(out_benchpdf, format='pdf', bbox_inches='tight')
    plt.close()
    df = bm.get_results(min_max_scale=False)
    print(df)
    df_transposed = df.transpose()
    df_transposed.to_csv(out_benchcsv, index=True)

if __name__ == "__main__":
    main()
