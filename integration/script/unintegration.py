### Date: 250812
### Image: harmony-py--
### Coder: ydgenomics

import matplotlib.pyplot as plt
import scanpy as sc
import pandas
import leidenalg
import click

@click.command()
@click.argument("input_h5ad", type=click.Path(exists=True))
@click.argument("out_h5ad", type=click.Path(exists=False), default=None)
@click.argument("out_umap", type=click.Path(exists=False), default=None)
@click.option('--batch_key', type=str, default=None, help="Batch key in identifying HVG and harmony integration")
@click.option('--sample_key', type=str, default=None, help="Sample key to distinguish sample")
@click.option('--cluster_key', type=str, default=None, help="Cluster key in samples")
@click.option('--resolution_set', type=float, default=None, help="set for resolution, is float")

def run_unintegration(input_h5ad, out_h5ad, out_umap, batch_key, sample_key, cluster_key, resolution_set):    
    # input_h5ad="/data/work/read/Cer_test.h5ad"
    orig_ad = sc.read_h5ad(input_h5ad)
    sc.pp.normalize_total(orig_ad, target_sum=1e4)
    sc.pp.log1p(orig_ad)
    #sc.pp.highly_variable_genes(orig_ad, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.highly_variable_genes(orig_ad)
    sc.pp.scale(orig_ad)
    sc.tl.pca(orig_ad, svd_solver="arpack")
    
    sc.pp.neighbors(orig_ad, n_neighbors=20, n_pcs=40)
    click.echo("computer cluster use leiden, and save the account of clusters in celltype")
    # resolution_set = 1.0
    sc.tl.leiden(orig_ad, resolution=resolution_set, key_added='celltype', flavor="igraph", n_iterations=2)
    sc.tl.umap(orig_ad)
    #batch_key="biosample"
    #sample_key="sample"
    #cluster_key="celltype"
    #out_umap="/data/work/unintegration/Cer_test_unintegration_integrated_UMAP.pdf"
    #out_h5ad="/data/work/unintegration/Cer_test_unintegration_integrated.h5ad"
    sc.pl.umap(orig_ad, color=[batch_key, cluster_key], legend_loc='on data', ncols=1)
    plt.savefig(out_umap, dpi=300,  bbox_inches='tight') 
    orig_ad.write(filename=out_h5ad,compression="gzip")

if __name__ == '__main__':
    run_unintegration()