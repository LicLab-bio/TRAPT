import json
import sys
from concurrent.futures import ThreadPoolExecutor

import anndata as ad
import numpy as np
import pandas as pd
from keras import Input, Model
from keras.layers import Dense
from keras.models import Model
from keras.optimizers import Adam
from numba import jit
from pyecharts import options as opts
from pyecharts.charts import Bar
from scipy import sparse
from tqdm import tqdm

from args import args_

args = args_(sys.argv[1])

@jit(nopython=True,nogil=True)
def get_f1_(params):
    i,j,y_true,y_pred = params
    beta=1
    tp = (y_true * y_pred).sum()
    fp = ((1 - y_true) * y_pred).sum()
    fn = (y_true * (1 - y_pred)).sum()
    epsilon = 1e-7
    precision = tp / (tp + fp + epsilon)
    recall = tp / (tp + fn + epsilon)
    f1 = (1 + beta**2)* (precision*recall) / (beta**2 * precision + recall + epsilon)
    return i,j,f1

@jit(nopython=True,nogil=True)
def get_auc_(params):
    i,j,labels,pred_prob,limit = params
    if np.any(labels) and np.sum(labels) > limit:
        index = np.argsort(pred_prob)
        l_p = labels > .5
        p_index = np.where(l_p)[0]
        n_index = np.where(~l_p)[0]
        l_p_rank = np.where(l_p[index])[0]
        rank_sum = np.sum(l_p_rank)
        auc = (rank_sum - len(p_index)*(1+len(p_index))/2) / (len(p_index)*len(n_index))
    else:
        auc = 0
    return i,j,auc

matrix = ad.read_h5ad('input/regu_matrix_tr_ad.h5ad')
data_ad = ad.read_h5ad('input/regu_matrix_atac_ad.h5ad')
samples = pd.read_csv('%s/sample_names_atac.csv' % args.output,header=None)[0].values
try:
    data = data_ad[samples].to_df().T.values
except:
    samples = list(map(lambda x:f'{x}_hg38',samples))
    data = data_ad[samples].to_df().T.values
input = Input(data.shape[1:])
x = Dense(args.dim, activation='sigmoid')(input)
output = Dense(1,activation='sigmoid')(x)
model = Model(input,output)
model.load_weights('%s/weight_atac.h5' % args.output)
gene_vec = model.predict(data,verbose=1).T[0]

def save_hist(data,bins,title):
    hist,bin_edges = np.histogram(data,bins=bins)
    bar=(
        Bar()
        .add_xaxis([str(x) for x in bin_edges[:-1]])
        .add_yaxis("",[float(x) for x in hist],category_gap=0,label_opts=opts.LabelOpts(is_show=False),)
        .set_global_opts(
            title_opts=opts.TitleOpts(title=title,pos_left="center"),
            legend_opts=opts.LegendOpts(is_show=False)
        )
    )
    bar.render(f"{args.output}/{title}.html")

bins = 100
hist,bin_edges = np.histogram(gene_vec,bins=bins)
save_hist(gene_vec,bins,'atac_rp_atac_pred')

auc_matrix = sparse.lil_matrix((matrix.shape[0],1))

def iter_params(matrix,gene_vec,trunk_size,trunk):
    start = trunk * trunk_size
    end = (trunk + 1) * trunk_size
    tr_dhs = matrix[start:end].to_df().values
    for i in range(tr_dhs.shape[0]):
        tr_vec = tr_dhs[i]
        yield start+i,0,gene_vec,tr_vec,args.min_inter_limit

with ThreadPoolExecutor(20) as pool:
    trunk_size = 512
    trunk_count = int(matrix.shape[0] / trunk_size) + 1
    for trunk in tqdm(range(trunk_count)):
        data = iter_params(matrix,gene_vec,trunk_size,trunk)
        tasks = [pool.submit(get_auc_, params) for params in data]
        for task in tasks:
            i,j,score = task.result()
            auc_matrix[i,j] = score

data = []
trs_base = list(matrix.obs.index.drop_duplicates())
for tr_base in tqdm(trs_base):
    trs_index = matrix.obs.loc[tr_base,'index']
    tr_samples = matrix.obs.loc[tr_base,'tr']
    if not bool(trs_index.shape):
        trs_index = [int(trs_index)]
        tr_samples = np.array([tr_samples])
    else:
        trs_index = list(trs_index)
        tr_samples = np.array(tr_samples)
    auc_vecs = auc_matrix[trs_index,:].toarray().flatten()
    rank_index = np.argsort(-auc_vecs)
    filter_top = rank_index[:args.top_n]
    tr_samples = tr_samples[filter_top]
    auc_vecs = auc_vecs[filter_top]
    data.append(auc_vecs.mean(axis=0))

data = np.vstack(data)
auc_matrix = sparse.csr_matrix(data)
auc_matrix_ad = ad.AnnData(data,obs={'tr':trs_base},dtype='float32')
auc_matrix_ad.write_h5ad('%s/auc_rp_atac_ad.h5ad' % args.output)
