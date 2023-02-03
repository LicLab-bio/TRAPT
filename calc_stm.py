import json
import os
import sys

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.nn.modules.loss
from scipy import sparse
from sklearn.preprocessing import normalize
from torch import optim
from torch.optim.lr_scheduler import StepLR

from args import args_

args = args_(sys.argv[1])

os.environ['CUDA_LAUNCH_BLOCKING'] = '1'

class GraphConvSparse(nn.Module):
    def __init__(self, input_dim, output_dim, adj, activation = F.relu, **kwargs):
        super(GraphConvSparse, self).__init__(**kwargs)
        self.weight = glorot_init(input_dim, output_dim)
        self.adj = adj
        self.activation = activation

    def forward(self, inputs):
        x = inputs
        x = torch.mm(x, self.weight)
        x = torch.mm(self.adj, x)
        outputs = self.activation(x)            
        return outputs

def sparse_to_tuple(sparse_mx):
    if not sparse.isspmatrix_coo(sparse_mx):
        sparse_mx = sparse_mx.tocoo()
    coords = np.vstack((sparse_mx.row, sparse_mx.col)).transpose()
    values = sparse_mx.data
    shape = sparse_mx.shape
    return coords, values, shape

def preprocess_graph(adj):
    adj = sparse.coo_matrix(adj)
    adj_ = adj + sparse.eye(adj.shape[0])
    rowsum = np.array(adj_.sum(1))
    degree_mat_inv_sqrt = sparse.diags(np.power(rowsum, -0.5).flatten())
    adj_normalized = adj_.dot(degree_mat_inv_sqrt).transpose().dot(degree_mat_inv_sqrt).tocoo()
    return adj_normalized

def glorot_init(input_dim, output_dim):
    init_range = np.sqrt(6.0/(input_dim + output_dim))
    initial = torch.rand(input_dim, output_dim)*2*init_range - init_range
    return nn.Parameter(initial)

class VGAE(nn.Module):
    def __init__(self, adj,input_dim,hidden1_dim,hidden2_dim):
        super(VGAE,self).__init__()
        self.base_gcn = GraphConvSparse(input_dim, hidden1_dim, adj)
        self.gcn_mean = GraphConvSparse(hidden1_dim, hidden2_dim, adj, activation=lambda x:x)
        self.gcn_stddev = GraphConvSparse(hidden1_dim, hidden2_dim, adj, activation=lambda x:x)

    def encode(self, X):
        hidden = self.base_gcn(X)
        self.mean = self.gcn_mean(hidden)
        self.std = self.gcn_stddev(hidden)
        self.mean[self.mean < 1e-32] = 1e-32
        self.std[self.std < 1e-32] = 1e-32
        if self.training:
            gaussian_noise = torch.randn_like(self.std)
            sampled_z = gaussian_noise.mul(self.std).add_(self.mean)
        else:
            sampled_z = self.mean
        return sampled_z
    def forward(self, X):
        Z = self.encode(X)
        A_pred = dot_product_decode(Z)
        return A_pred

def dot_product_decode(Z):
    A_pred = torch.matmul(Z,Z.t())
    A_pred = torch.sigmoid(A_pred)
    return A_pred

auc_tr_atac_ad = ad.read_h5ad('%s/auc_tr_atac_ad.h5ad' % args.output)
auc_tr_atac = auc_tr_atac_ad.X

auc_tr_chip_ad = ad.read_h5ad('%s/auc_tr_chip_ad.h5ad' % args.output)
auc_tr_chip = auc_tr_chip_ad.X

auc_tr_cre_ad = ad.read_h5ad('%s/auc_tr_cre_ad.h5ad' % args.output)
auc_tr_cre = auc_tr_cre_ad.X

auc_rp_atac_ad = ad.read_h5ad('%s/auc_rp_atac_ad.h5ad' % args.output)
auc_rp_atac = auc_rp_atac_ad.X

auc_rp_chip_ad = ad.read_h5ad('%s/auc_rp_chip_ad.h5ad' % args.output)
auc_rp_chip = auc_rp_chip_ad.X

with open(f'{args.output}/tr_top_samples.json','r') as file:
    samples_dict = json.loads(file.read())

samples = pd.DataFrame(samples_dict).T[0].values
data = pd.DataFrame(np.hstack([auc_tr_cre,auc_tr_atac,auc_tr_chip,auc_rp_atac,auc_rp_chip]).astype('float32'))

features = sparse.csr_matrix(data)

graph = ad.read_h5ad('input/graph.h5ad')
adj = graph.X

device = "cuda" if torch.cuda.is_available() and args.get_wc() == 0 else "cpu"
print(f"Using {device} device")
hidden1_dim = 32
hidden2_dim = 16
num_epoch = 500
learning_rate = 0.01
node_dim = features.shape[0]
input_dim = features.shape[1]
norm = node_dim * node_dim / float((node_dim * node_dim - adj.sum()) * 2)
adj_orig = adj
adj_orig = adj_orig - sparse.dia_matrix((adj_orig.diagonal()[np.newaxis, :], [0]), shape=adj_orig.shape)
adj_orig.eliminate_zeros()
adj_norm = preprocess_graph(adj)
adj_label = adj + sparse.eye(adj.shape[0])
adj_norm = sparse_to_tuple(adj_norm)
adj_label = sparse_to_tuple(adj_label)
features = sparse_to_tuple(features)
adj_norm = torch.sparse.FloatTensor(torch.LongTensor(adj_norm[0].T), 
                        torch.FloatTensor(adj_norm[1]), 
                        torch.Size(adj_norm[2])).to(device)
adj_label = torch.sparse.FloatTensor(torch.LongTensor(adj_label[0].T), 
                            torch.FloatTensor(adj_label[1]), 
                            torch.Size(adj_label[2])).to(device)
features = torch.sparse.FloatTensor(torch.LongTensor(features[0].T), 
                            torch.FloatTensor(features[1]), 
                            torch.Size(features[2])).to(device)

pos_weight = float(adj.shape[0] * adj.shape[0] - adj.sum()) / adj.sum()
weight_mask = adj_label.to_dense().view(-1) == 1
weight_tensor = torch.ones(weight_mask.size(0)) 
weight_tensor[weight_mask] = pos_weight
weight_tensor = weight_tensor.to(device)

STM = np.zeros(adj_orig.shape)
for i in range(args.vgae_iterations):
    print(f'# Model {i}')
    model = VGAE(adj_norm,input_dim,hidden1_dim,hidden2_dim).to(device)
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    scheduler = StepLR(optimizer, step_size=int(num_epoch/10), gamma=.99)
    print('# Train')
    for epoch in range(num_epoch):
        adj_pred = model(features)
        optimizer.zero_grad()
        logp = F.binary_cross_entropy_with_logits(adj_pred.view(-1), adj_label.to_dense().view(-1),weight=weight_tensor)
        p = torch.exp(-logp)
        loss = norm * (1 - p) ** 2 * logp
        kl_divergence = 0.5/node_dim * (1 + 2*torch.log(model.std).nan_to_num(0) - model.mean**2 - model.std**2).sum(1).mean()
        loss -= kl_divergence
        loss.backward()
        optimizer.step()
        scheduler.step()
        with torch.no_grad():
            if epoch % int(num_epoch/3) == 0:
                print("Epoch:", '%04d' % (epoch + 1), "train_loss=", "{:.5f}".format(loss.item()),'lr=',scheduler.get_last_lr())

    print('# Network reconfiguration')
    adj_pred = model(features)
    if device == 'cuda':
        STM += adj_pred.cpu().detach().numpy()
    else:
        STM += adj_pred.detach().numpy()

STM = STM/args.vgae_iterations

STM = sparse.csr_matrix(STM)

STM = STM - sparse.dia_matrix((STM.diagonal()[np.newaxis, :], [0]), shape=STM.shape)
STM = STM + sparse.eye(STM.shape[0])

from sklearn.metrics import auc, precision_recall_curve, roc_curve


def get_auc_auprc(label,pred):
    precision, recall, thresholds = precision_recall_curve(label, pred)
    auprc = auc(recall, precision)
    fpr, tpr, thersholds = roc_curve(label, pred)
    roc_auc = auc(fpr, tpr)
    return roc_auc,auprc

label = adj_orig.toarray().flatten() > .5
pred = STM.toarray().flatten()
roc_auc,auprc = get_auc_auprc(label,pred)
print(roc_auc,auprc)

sparse.save_npz('%s/STM.npz' % args.output,STM)
