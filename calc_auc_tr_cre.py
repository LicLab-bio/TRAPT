import json
import sys
from concurrent.futures import ThreadPoolExecutor

import anndata as ad
import numpy as np
import pandas as pd
from numba import jit
from scipy import sparse
from tqdm import tqdm

from args import args_

args = args_(sys.argv[1])

tr_dhs_ad = ad.read_h5ad('input/tr_filter_dhs_ad.h5ad')

chip_dhs = sparse.load_npz('%s/chip_dhs.npz' % args.output).toarray()
atac_dhs = sparse.load_npz('%s/atac_dhs.npz' % args.output).toarray()

chip_dhs = chip_dhs > .5
atac_dhs = atac_dhs > .5

sample_dhs = np.hstack([chip_dhs,atac_dhs]).astype(float)

significance_index_chip = pd.read_csv('%s/significance_index_chip.csv' % args.output,header=None)[0].values
significance_index_atac = pd.read_csv('%s/significance_index_atac.csv' % args.output,header=None)[0].values

significance_index = np.append(significance_index_chip,significance_index_atac)
data = pd.DataFrame(zip(significance_index,sample_dhs[0])).groupby(0).agg(max).reset_index()

sample_dhs = np.expand_dims(data[1].values,0)

tr_dhs_ad = tr_dhs_ad[:,data[0]]

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

auc_matrix = sparse.lil_matrix((tr_dhs_ad.shape[0],sample_dhs.shape[0]))

def iter_params(tr_dhs_ad,sample_dhs,trunk_size,trunk):
    start = trunk * trunk_size
    end = (trunk + 1) * trunk_size
    tr_dhs = tr_dhs_ad[start:end].to_df().values
    for i in range(tr_dhs.shape[0]):
        tr_vec = tr_dhs[i]
        for j in range(sample_dhs.shape[0]):
            sample_vec = sample_dhs[j]
            yield start+i,j,tr_vec,sample_vec,args.min_inter_limit

with ThreadPoolExecutor(20) as pool:
    trunk_size = 512
    trunk_count = int(tr_dhs_ad.shape[0] / trunk_size) + 1
    for trunk in tqdm(range(trunk_count)):
        data = iter_params(tr_dhs_ad,sample_dhs,trunk_size,trunk)
        tasks = [pool.submit(get_auc_, params) for params in data]
        for task in tasks:
            i,j,score = task.result()
            auc_matrix[i,j] = score

data = []
tr_top_samples = []
trs_base = list(tr_dhs_ad.obs.index.drop_duplicates())
for tr_base in tqdm(trs_base):
    trs_index = tr_dhs_ad.obs.loc[tr_base,'index']
    tr_samples = tr_dhs_ad.obs.loc[tr_base,'tr']
    if not bool(trs_index.shape):
        trs_index = [int(trs_index)]
        tr_samples = np.array([tr_samples])
    else:
        trs_index = list(trs_index)
        tr_samples = np.array(tr_samples)
    max_auc_vec = auc_matrix[trs_index,:].toarray().max(axis=0)
    auc_vecs = auc_matrix[trs_index,:].toarray().flatten()
    rank_index = np.argsort(auc_vecs)
    filter_top = rank_index[-args.top_n:]
    tr_samples = tr_samples[filter_top]
    tr_top_samples.append(list(tr_samples))
    auc_vecs = auc_vecs[filter_top]
    data.append(auc_vecs.mean(axis=0))

tr_top_samples = dict(zip(trs_base,tr_top_samples))
with open(f'{args.output}/tr_top_samples.json','w') as file:
    file.write(json.dumps(tr_top_samples))
    
data = np.vstack(data)
auc_matrix = sparse.csr_matrix(data)
auc_matrix_ad = ad.AnnData(data,obs={'tr':trs_base},dtype='float32')
auc_matrix_ad.write_h5ad('%s/auc_tr_cre_ad.h5ad' % args.output)
