import os
import sys

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse
from sklearn.preprocessing import normalize, scale

from args import args_

args = args_(sys.argv[1])

auc_tr_cre_ad = ad.read_h5ad('%s/auc_tr_cre_ad.h5ad' % args.output)
auc_tr_cre = auc_tr_cre_ad.to_df().values

auc_tr_atac_ad = ad.read_h5ad('%s/auc_tr_atac_ad.h5ad' % args.output)
auc_tr_atac = auc_tr_atac_ad.to_df().values

auc_tr_chip_ad = ad.read_h5ad('%s/auc_tr_chip_ad.h5ad' % args.output)
auc_tr_chip = auc_tr_chip_ad.to_df().values

auc_rp_atac_ad = ad.read_h5ad('%s/auc_rp_atac_ad.h5ad' % args.output)
auc_rp_atac = auc_rp_atac_ad.to_df().values

auc_rp_chip_ad = ad.read_h5ad('%s/auc_rp_chip_ad.h5ad' % args.output)
auc_rp_chip = auc_rp_chip_ad.to_df().values

tr = auc_tr_cre_ad.obs['tr'].values

batch = pd.DataFrame([1] * 3 + [2] * 2)
batch.to_csv('%s/batch.csv' % args.output,index=False,header=False)
data = pd.DataFrame(np.hstack([auc_tr_cre,auc_tr_atac,auc_tr_chip,auc_rp_atac,auc_rp_chip]).astype('float32'))
data.to_csv('%s/data.csv' % args.output,index=False,header=False)
os.system('Rscript sva.r %s' % args.output)
combat_data = pd.read_csv('%s/combat_data.csv' % args.output).values

STM = sparse.load_npz('%s/STM.npz' % args.output).toarray()
interaction_score = np.expand_dims(np.mean(STM,axis=1),1)
auc_score = np.expand_dims(np.mean(combat_data,axis=1),1)
tr_activity = np.multiply(interaction_score,auc_score)

activity_summary = pd.DataFrame(np.hstack([tr_activity,interaction_score,data]),index=tr,columns=['tr_activity','interaction_score','re_score','atac_score','chip_score','rp_atac_score','rp_chip_score'])
activity_summary.sort_values('tr_activity',ascending=False).to_csv('%s/activity_summary.csv' % args.output) 
