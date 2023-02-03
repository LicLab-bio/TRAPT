import os
import subprocess
import sys

import anndata
import numpy as np
import pandas as pd
from keras import Input, Model
from keras import backend as K
from keras.layers import Dense, Dropout
from keras.models import Model
from keras.optimizers import Adam
from sklearn import metrics

from args import args_

args = args_(sys.argv[1])
sample = sys.argv[2]

if args.get_wc() > 0:
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    os.environ["CUDA_VISIBLE_DEVICES"] = ""

def get_auc_auprc(label,pred):
    precision, recall, _ = metrics.precision_recall_curve(label, pred)
    auprc = metrics.auc(recall, precision)
    fpr, tpr, _ = metrics.roc_curve(label, pred)
    roc_auc = metrics.auc(fpr, tpr)
    return roc_auc,auprc

data_ad = anndata.read_h5ad(f'input/regu_matrix_{sample}_ad.h5ad')
geneset = pd.read_csv(args.input,header=None)[0]
genes = data_ad.var_names.values
sample_names = data_ad.obs_names.values
status = np.zeros(len(genes))
mask = np.in1d(genes,geneset)

pos_weight = float(len(mask) - mask.sum()) / mask.sum()
weight_mask = mask == 1
sample_weight = np.ones(len(weight_mask))
sample_weight[weight_mask] = pos_weight

X = data_ad.T.to_df().values
T = mask.astype(int)

lamba = {
    "chip":0.005,
    "atac":0.005
}

# Sparsity constraint
L1 = lambda weight_matrix:lamba[sample] * K.sum(K.sqrt(K.tf.reduce_sum(K.square(weight_matrix), axis=1)))

# The teacher network input (X,T) is used to generate the low-dimensional representation Y.
input = Input(X.shape[1:])
input = Dropout(.2)(input)
y = Dense(args.dim, activation='sigmoid')(input)
feature_extraction = Model(input,y)
t = feature_extraction(input)
t = Dropout(.2)(t)
output = Dense(1,activation='sigmoid')(t)
teacher = Model(input,output)
teacher.compile(optimizer=Adam(learning_rate=lamba[sample]/10), loss='mse')
teacher.summary()
teacher.fit(X, T, epochs=64, batch_size=256, sample_weight=sample_weight,verbose=0)

T_pred = teacher.predict(X)
print('Teacher AUC,AUPRC is ',get_auc_auprc(T, T_pred))

Y = feature_extraction.predict(X)

# Students input (X,Y) into the network to get the first layer weight ranking.
input = Input((X.shape[-1],))
x = Dense(Y.shape[-1] * 4, activation='sigmoid', kernel_regularizer=L1)(input)
x = Dropout(.2)(x)
y = Dense(Y.shape[-1], activation='sigmoid')(x)
student = Model(input,y)
student.compile(optimizer=Adam(learning_rate=lamba[sample]/10), loss='mse')
student.summary()
student.fit(X, Y, epochs=64, batch_size=256,sample_weight=sample_weight,verbose=0)

# get the first layer weights.
weights = student.layers[1].get_weights()[0]
# get the feature importance.
feature_importances = np.sum(np.square(weights),1)

feature_indexs = np.argsort(-feature_importances)

for i in range(len(feature_indexs)):
    indexs = feature_indexs[:i]
    if i > 50:
        feature_indexs = feature_indexs[:10]
        print("feature_indexs",feature_indexs.shape)
        break
    if feature_importances[indexs].sum() > .99 * feature_importances.sum():
        feature_indexs = indexs
        print("feature_indexs",feature_indexs.shape)
        break

X_filter = X[:,feature_indexs]
# The new model is used to predict the enhanced sub-signal.
input = Input(X_filter.shape[1:])
input = Dropout(.2)(input)
x = Dense(args.dim, activation='sigmoid')(input)
x = Dropout(.2)(x)
output = Dense(1,activation='sigmoid')(x)
model = Model(input,output)
model.compile(optimizer=Adam(learning_rate=lamba[sample]/10), loss='mse')
model.summary()
model.fit(X_filter, T, epochs=64, batch_size=256,sample_weight=sample_weight,verbose=0)

T_pred = model.predict(X_filter)
print('Filter AUC,AUPRC is ',get_auc_auprc(T, T_pred))

model.save_weights(filepath=f'{args.output}/weight_{sample}.h5')
pd.DataFrame(sample_names[feature_indexs]).to_csv(f'{args.output}/sample_names_{sample}.csv',index=False,header=False)
