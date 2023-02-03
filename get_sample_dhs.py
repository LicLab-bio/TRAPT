import os
import sys

import numpy as np
import pandas as pd
from keras import Input, Model
from keras.layers import Dense, Dropout
from keras.models import Model
from pyecharts import options as opts
from pyecharts.charts import Bar
from scipy import sparse
from sklearn.metrics import auc, precision_recall_curve, roc_curve

from args import args_

args = args_(sys.argv[1])
sample = sys.argv[2]

if args.get_wc() > 0:
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    os.environ["CUDA_VISIBLE_DEVICES"] = ""

def get_auc_auprc(label,pred):
    precision, recall, thresholds = precision_recall_curve(label, pred)
    auprc = auc(recall, precision)
    fpr, tpr, thersholds = roc_curve(label, pred)
    roc_auc = auc(fpr, tpr)
    return roc_auc,auprc

samples = pd.read_csv(f'{args.output}/sample_names_{sample}.csv',header=None)[0]
sample_dhs = []
for sample_ in samples:
    print(sample_)
    vec = pd.read_csv(f'input/{sample}_matrix/dhs@{sample_}.csv',sep='\t',header=None)
    vec = vec.iloc[:,0].values
    sample_dhs.append(vec)

sample_dhs = np.vstack(sample_dhs).T

input = Input(sample_dhs.shape[1:])
x = Dense(args.dim, activation='sigmoid')(input)
output = Dense(1,activation='sigmoid')(x)
model = Model(input,output)
model.load_weights(f'{args.output}/weight_{sample}.h5')

sample_dhs_pred = model.predict(sample_dhs,verbose=0)

def save_hist(sample_dhs_filter,bins,title):
    hist,bin_edges = np.histogram(sample_dhs_filter,bins=bins)
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

bins = 1000
input_gene_dhs_index = args.get_dhs_index()
hist,bin_edges = np.histogram(sample_dhs_pred[input_gene_dhs_index],bins=bins)
save_hist(sample_dhs_pred[input_gene_dhs_index],bins,f'{sample}_dhs_pred')

# evs = bin_edges[signal.argrelextrema(hist, np.greater)]

for i in range(1,len(hist)):
    if hist[:i].sum() > args.range_threshold * hist.sum():
        down = bin_edges[i]
        print('# # # down # # #',down)
        break
for i in range(1,len(hist)):
    if hist[-i:].sum() > args.range_threshold * hist.sum():
        up = bin_edges[-i]
        print('# # # up # # #',up)
        break
significance_index = np.where(np.any(np.hstack([
    sample_dhs_pred[input_gene_dhs_index] <= down,
    sample_dhs_pred[input_gene_dhs_index] >= up
]),axis=1))[0]
significance_index = input_gene_dhs_index[significance_index]
sample_dhs_filter = sample_dhs_pred[significance_index]
save_hist(sample_dhs_filter,bins,f'{sample}_dhs_filter')

pd.DataFrame(significance_index).to_csv(f'{args.output}/significance_index_{sample}.csv',index=False,header=False)
sparse.save_npz(f'{args.output}/{sample}_dhs.npz',sparse.csr_matrix(sample_dhs_filter).T)
sparse.save_npz(f'{args.output}/{sample}_dhs_pred.npz',sparse.csr_matrix(sample_dhs_pred).T)
