# dealplus.python
import click
import matplotlib.pyplot as plt
import scanpy as sc
from matplotlib.backends.backend_pdf import PdfPages

@click.command()
@click.argument("input_h5ad", type=click.Path(exists=True))
@click.argument("out_h5ad", type=click.Path(exists=False), default=None)
@click.argument("out_umap", type=click.Path(exists=False), default=None)
@click.option('--other1_key', type=str, default=None, help="Sample key to distinguish sample")
@click.option('--other2_key', type=str, default=None, help="Cluster key in samples")

def run_dealplus(input_h5ad, out_h5ad, out_umap, other1_key, other2_key):
    sc.set_figure_params(dpi_save=300, frameon=False, figsize=(10, 6))
    adata = sc.read_h5ad(input_h5ad)
    #add biosample.celltype column
    adata.obs['biosample.celltype'] = adata.obs['biosample'].astype(str) + '_' + adata.obs['celltype'].astype(str)
    print(adata.obs["biosample.celltype"].unique())
    #visual
    with PdfPages(out_umap) as pdf:
        sc.pl.umap(adata, color=["biosample", "sample", "celltype", "biosample.celltype"], ncols=1)
        plt.savefig(pdf, format='pdf', dpi=300, bbox_inches='tight')
        plt.close()
        required_cols = ['total_counts', 'n_genes']
        if all(col in adata.obs.columns for col in required_cols):
            sc.pl.violin(adata, keys=['total_counts'], log=False, groupby="biosample", show=False)
            plt.savefig(pdf, format='pdf', dpi=300, bbox_inches='tight')
            plt.close()
            sc.pl.violin(adata, keys=['n_genes'], log=False, groupby="biosample", show=False)
            plt.savefig(pdf, format='pdf', dpi=300, bbox_inches='tight')
            plt.close()
        else:
            print("obs lacked total_counts or n_genes column")
    # 检查 obs 中是否存在这两个列
    if other1_key in adata.obs.columns and other2_key in adata.obs.columns:
        # 检查是否满足条件：other1_key 不等于 "biosample" 且 other2_key 不等于 "celltype"
        if other1_key != "biosample" or other2_key != "celltype":
            # 添加组合列
            adata.obs['other1_other2'] = adata.obs[other1_key].astype(str) + '_' + adata.obs[other2_key].astype(str)

            # 打印组合列的唯一值
            print(adata.obs["other1_other2"].unique())

            # 绘制 UMAP 并保存为 PDF
            with PdfPages('other'+out_umap) as pdf:
                sc.pl.umap(adata, color=[other1_key, "sample", other2_key, "other1_other2"], ncols=1, show=False)
                plt.savefig(pdf, format='pdf', dpi=300, bbox_inches='tight')
                plt.close()
        else:
            print("other1_key 或 other2_key 为默认值，跳过 UMAP 绘图。")
    else:
        print("obs 中缺少指定的列，跳过 UMAP 绘图。")
    
    click.echo("Save output")
    adata.write(filename=out_h5ad,compression="gzip")


if __name__ == '__main__':
    run_dealplus()
