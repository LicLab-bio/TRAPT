import functools
import os
import sys

import numpy as np
import pandas as pd

from TRAPT.Tools import RP_Matrix, Args, Type
from TRAPT.CalcTRAUC import CalcTRAUC
from TRAPT.DLFS import FeatureSelection


def get_params(func):
    @functools.wraps(func)
    def wrapper(args):
        return func(*args)

    return wrapper


@get_params
def runTRAPT(rp_matrix: RP_Matrix, args: Args):
    obs = rp_matrix.TR.obs

    if os.path.exists(f"{args.output}/H3K27ac_RP.csv"):
        H3K27ac_RP = pd.read_csv(f"{args.output}/H3K27ac_RP.csv", header=None)[0]
    else:
        FS_H3K27ac = FeatureSelection(args, rp_matrix.H3K27ac, Type.H3K27ac)
        H3K27ac_RP = FS_H3K27ac.run()
        H3K27ac_RP.to_csv(f"{args.output}/H3K27ac_RP.csv", index=False, header=False)

    if os.path.exists(f"{args.output}/ATAC_RP.csv"):
        ATAC_RP = pd.read_csv(f"{args.output}/ATAC_RP.csv", header=None)[0]
    else:
        FS_ATAC = FeatureSelection(args, rp_matrix.ATAC, Type.ATAC)
        ATAC_RP = FS_ATAC.run()
        ATAC_RP.to_csv(f"{args.output}/ATAC_RP.csv", index=False, header=False)

    if os.path.exists(f"{args.output}/RP_TR_H3K27ac_auc.csv"):
        RP_TR_H3K27ac_auc = pd.read_csv(
            f"{args.output}/RP_TR_H3K27ac_auc.csv", index_col=0, header=None
        )
    else:
        H3K27ac_RP = H3K27ac_RP.values.flatten()
        CTR_TR = CalcTRAUC(args, rp_matrix.TR_H3K27ac, H3K27ac_RP)
        RP_TR_H3K27ac_auc = CTR_TR.run()
        RP_TR_H3K27ac_auc.to_csv(f"{args.output}/RP_TR_H3K27ac_auc.csv", header=False)

    if os.path.exists(f"{args.output}/RP_TR_ATAC_auc.csv"):
        RP_TR_ATAC_auc = pd.read_csv(
            f"{args.output}/RP_TR_ATAC_auc.csv", index_col=0, header=None
        )
    else:
        ATAC_RP = ATAC_RP.values.flatten()
        CTR_TR = CalcTRAUC(args, rp_matrix.TR_ATAC, ATAC_RP)
        RP_TR_ATAC_auc = CTR_TR.run()
        RP_TR_ATAC_auc.to_csv(f"{args.output}/RP_TR_ATAC_auc.csv", header=False)

    data_auc = pd.concat([RP_TR_H3K27ac_auc, RP_TR_ATAC_auc], axis=1)
    data_auc /= np.linalg.norm(data_auc, axis=0, keepdims=True)
    TR_activity = pd.DataFrame(
        np.sum(data_auc.values, axis=1), index=data_auc.index, columns=[1]
    )
    TR_detail = pd.concat([TR_activity, data_auc], axis=1).reset_index()
    TR_detail.columns = ["TR", "TR activity", "RP_TR_H3K27ac_auc", "RP_TR_ATAC_auc"]
    obs.index.name = "TR"
    TR_detail = TR_detail.merge(obs.reset_index(), on="TR").sort_values(
        "TR activity", ascending=False
    )
    TR_detail.to_csv(os.path.join(args.output, "TR_detail.txt"), index=False, sep="\t")
    return TR_detail


if __name__ == "__main__":
    input = sys.argv[1]
    output = sys.argv[2]
    library = sys.argv[3]
    rp_matrix = RP_Matrix(library)
    args = Args(input, output)
    os.system(f"mkdir -p {output}")
    runTRAPT([rp_matrix, args])
