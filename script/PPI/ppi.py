import os
from optparse import OptionParser

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse

usage = "usage: %prog [options] -p [PATH] -o [OUTPUT]"
parser = OptionParser(usage = usage)
parser.add_option("-p", dest="PATH",nargs = 1, default=None,help = "TCAPT tool path.")
parser.add_option("-o", dest="OUTPUT",nargs = 1, default=None,help = "Enter an output folder.")
options,args = parser.parse_args()

def draw_heatmap(matrix, name):
    if isinstance(matrix,sparse.data._data._data_matrix):
        matrix_ = matrix.toarray()
    elif isinstance(matrix,list):
        matrix_ = np.asanyarray(matrix)
    else:
        matrix_ = matrix
    plt.figure(figsize=(60,50))
    sns.heatmap(pd.DataFrame(matrix_))
    plt.savefig(name)
    plt.clf()

def table_to_sparse(ind_dict,table,shape):
    table = table.iloc[:,:2]
    table.columns = ['row_ind','col_ind']
    table = table[table.apply(lambda x:x['row_ind'] in ind_dict and x['col_ind'] in ind_dict,axis=1)]
    table['row_ind'] = table['row_ind'].apply(ind_dict.get)
    table['col_ind'] = table['col_ind'].apply(ind_dict.get)
    print(table)
    add_row = table[['col_ind','row_ind']]
    add_row.columns = ['row_ind','col_ind']
    table = pd.concat([table,add_row],ignore_index=True)
    table = table.drop_duplicates()
    table['data'] = 1
    print(table)
    table = table.sort_values(['col_ind'],ascending=True,ignore_index=True)
    matrix = sparse.csc_matrix((table['data'],(table['row_ind'],table['col_ind'])),shape=(shape,shape),dtype='float32')
    print(matrix.shape)
    return matrix

def get_graph(file,tr_name):
    table = pd.read_csv(file)
    TR_dict = dict(zip(tr_name,range(len(tr_name))))
    matrix = table_to_sparse(TR_dict,table,len(tr_name))
    graph = ad.AnnData(X=matrix)
    graph.obs_names = tr_name
    graph.var_names = tr_name
    return graph

tr_dhs_ad = ad.read_h5ad(f'{options.PATH}/input/tr_filter_dhs_ad.h5ad')
tr_dhs_ad.obs = tr_dhs_ad.obs.reset_index()
tr_dhs_ad.obs['tr_base'] = tr_dhs_ad.obs['tr'].apply(lambda x:x.split('@')[0])
tr_dhs_ad.obs = tr_dhs_ad.obs.set_index('tr_base')
tr_dhs_ad.obs['index'] = tr_dhs_ad.obs['index'].apply(int)
trs_base = list(tr_dhs_ad.obs.index.drop_duplicates())

graph = get_graph('annotation/human_ppi.csv',trs_base)
if not os.path.exists(options.OUTPUT):
    os.mkdir(options.OUTPUT)

graph.write_h5ad(f'{options.OUTPUT}/graph.h5ad')

draw_heatmap(graph.to_df(), f'{options.OUTPUT}/graph.png')
