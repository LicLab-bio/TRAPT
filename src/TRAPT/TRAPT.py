import functools
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler, StandardScaler

from TRAPT.CalcTRAUC import CalcTRAUC
from TRAPT.DLFS import FeatureSelection


def get_params(func):
    @functools.wraps(func)
    def wrapper(args):
        return func(*args)

    return wrapper


class Args:
    def __init__(self, input, output, threads=16, trunk_size=2048 * 16) -> None:
        self.input = input
        self.output = output
        self.threads = threads
        self.trunk_size = trunk_size


class Type:
    H3K27ac = 'H3K27ac'
    ATAC = 'ATAC'


class RPMatrix:
    def __init__(self, library, name, to_array=True):
        self.data = ad.read_h5ad(os.path.join(library, name))
        if to_array:
            self.data.X = self.data.to_df().values

    def norm(self, type="l2", axis=1):
        assert type in ["l1", "l2"]
        self.data.X = self.data.X - self.data.X.min()
        if type == "l1":
            self.data.X /= np.sum(self.data.X, axis=axis, keepdims=True).clip(min=1e-17)
        if type == "l2":
            self.data.X /= np.linalg.norm(self.data.X, axis=axis, keepdims=True).clip(
                min=1e-17
            )
        return self

    def standard_scale(self, axis=1):
        ss = StandardScaler()
        if axis == 0:
            self.data.X = ss.fit_transform(self.data.X)
        if axis == 1:
            self.data.X = ss.fit_transform(self.data.X.T).T
        return self

    def binarization(self):
        self.data.X = (self.data.X > 1e-17).astype(np.float32)
        return self

    def minmax_scale(self, axis=1):
        ss = MinMaxScaler()
        if axis == 0:
            self.data.X = ss.fit_transform(self.data.X)
        if axis == 1:
            self.data.X = ss.fit_transform(self.data.X.T).T
        return self

    def add(self, data):
        self.data.X += data.X
        return self

    def get_data(self):
        return self.data


class RP_Matrix:
    def __init__(self, library) -> None:
        self.TR = RPMatrix(library, 'RP_Matrix_TR.h5ad').norm().get_data()
        self.TR_H3K27ac = (
            RPMatrix(library, 'RP_Matrix_TR_H3K27ac.h5ad')
            .norm()
            .add(self.TR)
            .get_data()
        )
        self.TR_ATAC = (
            RPMatrix(library, 'RP_Matrix_TR_ATAC.h5ad').norm().add(self.TR).get_data()
        )
        self.H3K27ac = (
            RPMatrix(library, 'RP_Matrix_H3K27ac.h5ad').standard_scale().get_data()
        )
        self.ATAC = RPMatrix(library, 'RP_Matrix_ATAC.h5ad').standard_scale().get_data()


@get_params
def runTRAPT(rp_matrix: RP_Matrix, args: Args):
    obs = rp_matrix.TR.obs

    if os.path.exists(f'{args.output}/H3K27ac_RP.csv'):
        H3K27ac_RP = pd.read_csv(f'{args.output}/H3K27ac_RP.csv', header=None)[0]
    else:
        FS_H3K27ac = FeatureSelection(args, rp_matrix.H3K27ac, Type.H3K27ac)
        H3K27ac_RP = FS_H3K27ac.run()
        H3K27ac_RP.to_csv(f'{args.output}/H3K27ac_RP.csv', index=False, header=False)

    if os.path.exists(f'{args.output}/ATAC_RP.csv'):
        ATAC_RP = pd.read_csv(f'{args.output}/ATAC_RP.csv', header=None)[0]
    else:
        FS_ATAC = FeatureSelection(args, rp_matrix.ATAC, Type.ATAC)
        ATAC_RP = FS_ATAC.run()
        ATAC_RP.to_csv(f'{args.output}/ATAC_RP.csv', index=False, header=False)

    if os.path.exists(f'{args.output}/RP_TR_H3K27ac_auc.csv'):
        RP_TR_H3K27ac_auc = pd.read_csv(
            f'{args.output}/RP_TR_H3K27ac_auc.csv', index_col=0, header=None
        )
    else:
        H3K27ac_RP = H3K27ac_RP.values.flatten()
        CTR_TR = CalcTRAUC(args, rp_matrix.TR_H3K27ac, H3K27ac_RP)
        RP_TR_H3K27ac_auc = CTR_TR.run()
        RP_TR_H3K27ac_auc.to_csv(f'{args.output}/RP_TR_H3K27ac_auc.csv', header=False)

    if os.path.exists(f'{args.output}/RP_TR_ATAC_auc.csv'):
        RP_TR_ATAC_auc = pd.read_csv(
            f'{args.output}/RP_TR_ATAC_auc.csv', index_col=0, header=None
        )
    else:
        ATAC_RP = ATAC_RP.values.flatten()
        CTR_TR = CalcTRAUC(args, rp_matrix.TR_ATAC, ATAC_RP)
        RP_TR_ATAC_auc = CTR_TR.run()
        RP_TR_ATAC_auc.to_csv(f'{args.output}/RP_TR_ATAC_auc.csv', header=False)

    data_auc = pd.concat([RP_TR_H3K27ac_auc, RP_TR_ATAC_auc], axis=1)
    data_auc /= np.linalg.norm(data_auc, axis=0, keepdims=True)
    TR_activity = pd.DataFrame(
        np.sum(data_auc.values, axis=1), index=data_auc.index, columns=[1]
    )
    TR_detail = pd.concat([TR_activity, data_auc], axis=1).reset_index()
    TR_detail.columns = ['TR', 'TR activity', 'RP_TR_H3K27ac_auc', 'RP_TR_ATAC_auc']
    obs.index.name = 'TR'
    TR_detail = TR_detail.merge(obs.reset_index(), on='TR').sort_values(
        'TR activity', ascending=False
    )
    TR_detail.to_csv(os.path.join(args.output, 'TR_detail.txt'), index=False, sep='\t')
    return TR_detail


if __name__ == '__main__':
    input = sys.args[0]
    output = sys.args[1]
    library = sys.args[2]
    rp_matrix = RP_Matrix(library)
    args = Args(input, output)
    os.system(f'mkdir -p {output}')
    runTRAPT([rp_matrix, args])
