import plotly.graph_objects as go
import scanpy as sc
import pandas as pd
import click
import itertools
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
from operator import itemgetter

def get_tbs(prefix,sequences):
    n=0
    tbs = {}
    while n+1<len(sequences):
        path1 = prefix + '_' + sequences[n] + '_' + sequences[n+1] + '.h5ad'
        path2 = prefix + '_' + sequences[n+1] + '_' + sequences[n] + '.h5ad'
        if os.path.exists(path1):
            path=path1
        elif os.path.exists(path2):
            path=path2
        adata = sc.read_h5ad(path)
        temp=adata.obs[['celltype','logit_inferred','sccaf_projection']]
        spe=sequences[n] + '_from'
        temp1=temp[temp.sccaf_projection.str.contains(spe)]
        temp1.loc[:,'type']=temp1['celltype'].astype(str)+'_to_'+temp1['logit_inferred'].astype(str)
        temp2=temp1['type'].value_counts()
        temp3=temp2.index.map(lambda x:x.split('_to_')[0])
        temp4=temp2.index.map(lambda x:x.split('_to_')[1])
        data={'values':temp2,'type1':temp3,'type2':temp4}
        temp5=pd.DataFrame(data)
        tbs[n]=temp5
        n=n+1
    return tbs


def sanky_plot(tbs):
    temp=pd.concat(tbs,axis=0)
    types=pd.unique(np.concatenate((temp['type1'].values,temp['type2'].values)))
    indexs=dict(zip(types, range(len(types))))
    colors=plt.cm.plasma(np.linspace(0, 1, len(indexs)))
    fig = go.Figure(data=[go.Sankey(
    node = dict(
      pad = 15,
      thickness = 20,
      line = dict(color = "black", width = 0.5),
      label = types,
    ),
    link = dict(
      source = list(itemgetter(*temp['type1'].values)(indexs)), # indices correspond to labels, eg A1, A2, A1, B1, ...
      target = list(itemgetter(*temp['type2'].values)(indexs)),
      value = list(temp['values'].values)
    ))])
    return fig

@click.command()
@click.argument("sequence", type=click.Path(exists=True))
@click.option('--prefix', type=str, default=None, help="prefix of target h5ad files")

def sanky_plot_2(prefix,sequence):
    sequences=[]
    with open(sequence,'r') as f:
        for line in f:
             elem = ''.join(line.strip('\n').split(','))
             sequences.append(elem)


    tbs1=get_tbs(prefix,sequences)
    tbs2=get_tbs(prefix,sequences[::-1])
    fig1=sanky_plot(tbs1)
    fig2=sanky_plot(tbs2)
    fig1.write_html(prefix+"_1.html")
    fig2.write_html(prefix+"_2.html") 

if __name__ == '__main__':
    sanky_plot_2()
